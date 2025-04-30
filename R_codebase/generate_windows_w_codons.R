require(GenomicFeatures)

## ----- helper functions for filtering and annotating windows -----

gr_from_bed <- function(bedfile) {
  bed <-
    read.delim(
      bedfile,
      sep = "\t",
      header = F,
      stringsAsFactors = F
    )
  if(ncol(bed) == 6){
    colnames(bed) <-
      c('chr', 'start', 'end', 'id', 'score', 'strand')
  }
  if(ncol(bed) == 4){
    colnames(bed) <-
      c('chr', 'start', 'end', 'id')
  }
  
  bed.gr <-
    makeGRangesFromDataFrame(
      bed,
      keep.extra.columns = T,
      starts.in.df.are.0based = T,
      seqinfo = Seqinfo(seqnames = as.character(unique(bed$chr)))
    )
  return(bed.gr)
}

mygetSeq <- ## update seqlevels before calling getSeq
  function(tx.seqs, gr, ...){
  seqlevels(gr) <- seqlevelsInUse(gr)
  getSeq(tx.seqs, gr, ...)
}

win_filter <- ## filter wins that fall off transcript boundaries
  function(win, txmodel = txdat$txmodel){
  win <- win[start(win) >1] #filter windows before TSS
  win <- win[end(win) < txmodel[match(as.character(seqnames(win)), txmodel$tx_name),]$tx_len] # filter windows after TES
  return(win)
}

get_frame <- ## get the reading frame of the first nucleotide of a GRanges in transcriptome space
  function(gr, txmodel = txdat$txmodel){
    pos <- start(gr)   # beware of 1-based vs. 0-based
    ## need to filter sites outside of CDS, as the modulo approach will otherwise assign a frame
    pos[pos <  txmodel[as.integer(match(seqnames(gr), txmodel$tx_name)),]$cds.start]   <- NA
    pos[pos >= txmodel[as.integer(match(seqnames(gr), txmodel$tx_name)),]$utr3.start ] <- NA
    ## assign frame using modulo of position in CDS
    frame <-   (pos - txmodel[as.integer(match(seqnames(gr), txmodel$tx_name)),]$cds.start) %% 3 # beware of 1-based vs. 0-based
    return(factor(frame))
  }

## a function to convert internal frame (modulo 1st nt of start codon) to print version (position in codon)
get_framelabel <- function(frame){
  factor(ifelse(is.na(frame), as.character(NA), paste("Frame", as.numeric(as.character(frame))+1)))
}


get_codon <- ## assign the first nucleotide of a GRanges in transcriptome space with the codon it hits
  function(gr, txmodel = txdat$txmodel, txs = tx.seqs){
    pos <- start(gr)   # beware of 1-based vs. 0-based
    ## need to filter sites outside of CDS, as the modulo approach will otherwise assign a frame
    pos[pos <  txmodel[as.integer(match(seqnames(gr), txmodel$tx_name)),]$cds.start]   <- NA
    pos[pos >= txmodel[as.integer(match(seqnames(gr), txmodel$tx_name)),]$utr3.start ] <- NA
    ## assign frame using modulo of position in CDS
    frame <-   (pos - txmodel[as.integer(match(seqnames(gr), txmodel$tx_name)),]$cds.start) %% 3 # beware of 1-based vs. 0-based
    ## work only on gr that has an annotated frame
    has_frame <- !is.na(frame)
    in_txs <- seqnames(gr) %in% names(txs)
    # assign codon vector with NA to be in parallel with gr
    codon <- as.character(rep(NA,length(gr)))
    codon[as.logical(has_frame & in_txs)] <-
      as.character(mygetSeq(
       txs,
       resize(GenomicRanges::shift(gr[as.logical(has_frame & in_txs)], -frame[as.logical(has_frame & in_txs)]), 3 ) 
       #gr[as.logical(has_frame & in_txs)] 
      )
      )
    return(codon)
  }

get_codon_repeats <-
  function(gr = NULL,
           cdn = NULL,
           repmax = 5) {
    as.data.table(foreach(
      i = 1:repmax,
      .combine = cbind,
      .export = "get_kmer",
      .packages = c("GenomicRanges")
    ) %dopar% {
      ifelse(get_kmer(resize(gr, i * 3, "start")) == paste(rep(cdn, i), collapse =
                                                             ""), i, 0)
    })[, do.call(pmax, .SD)]
  }



get_codon_gr <-
  function(txmodel = NULL, tx.seqs = NULL, eejs.tx = NULL, drach.tx = NULL, ex_on_tx = NULL){
    
    cds.gr <-
      txmodel[,GRanges(seqnames = tx_name, ranges =  IRanges(start=cds.start, end=utr3.start-1))] ## generate GRanges from txmodel CDS annotation
    
    cds.gr <-
      cds.gr[(width(cds.gr) %% 3) == 0]  ## remove transcripts with CDS not divisible by 3
    
    ### generate codon GRanges
    cods.gr <- unlist(tile(cds.gr, width = 3))
    strand(cods.gr) <- "+"
    
    ### unique codonID (used for finding unmethylated control codons)
    cods.gr$codonID <-
      data.table(tx_name = as.character(seqnames(cods.gr)))[, codonID := paste0(gsub("\\.*$","",tx_name),"_", seq_len(.N)), by = "tx_name"]$codonID
    
    cods.gr$codon <-
      get_kmer(gr = cods.gr, txs = tx.seqs)
    
    cods.gr$kmer <-
      get_kmer(cods.gr+2, tx.seqs)
    
    ### codon position in CDS (relative)
    cods.gr$relpos.cds <-
      relpos.cds.gr(cods.gr)  
    
    ### codon position in CDS (absolute)
    cods.gr$codon_rank <-
      data.table(tx_name = as.character(seqnames(cods.gr)), pos = as.integer(start(cods.gr))  )[, rank := order(pos), by=c("tx_name")][, rank]
    
    cods.gr$is_final_exon <-
      (start(cods.gr) - txmodel$exex_pos_last[as.integer(match(seqnames(cods.gr), txmodel$tx_name  ))]) > 0 
    
    ### DRACH codons
    cods.gr$is_drach <-
      overlapsAny(cods.gr, drach.tx) ## beware that for ACA codons, this might result in the codon being called as DRACH in the unmethylatable frame2!
    #table(overlapsAny(resize(cods.gr[cods.gr$codon == "ACA"],1,"end"), drach.tx))
    
    ### define the frame in which a putative m6A might be
    cods.gr$m6A_frame <- rep(as.integer(NA), length(cods.gr))
    cods.gr[cods.gr$codon %in% c("ACT", "ACA", "ACC")]$m6A_frame = 0
    cods.gr[cods.gr$codon %in% c("GAC", "AAC")]$m6A_frame = 1
    cods.gr[cods.gr$codon %in% c("AAA", "AGA", "GAA", "GGA")]$m6A_frame = 2

    ### EEJ distance
    cods.gr <- resize(cods.gr, 1, "start")
    
    cods.gr$eejFollow <-  follow(cods.gr, eejs.tx)
    cods.gr$eejFollowDist <- as.integer(NA)
    cods.gr[!is.na(cods.gr$eejFollow)]$eejFollowDist <-
      start(cods.gr[!is.na(cods.gr$eejFollow)]) - start(eejs.tx[cods.gr[!is.na(cods.gr$eejFollow)]$eejFollow]) -1
    cods.gr$eejFollowName <- as.character(NA)
    cods.gr[!is.na(cods.gr$eejFollow)]$eejFollowName <-
      eejs.tx[cods.gr[!is.na(cods.gr$eejFollow)]$eejFollow]$exon_name
    
    cods.gr$eejPrecede <-  precede(cods.gr, eejs.tx)
    cods.gr$eejPrecedeDist <- as.integer(NA)
    cods.gr[!is.na(cods.gr$eejPrecede)]$eejPrecedeDist <-
      start(cods.gr[!is.na(cods.gr$eejPrecede)]) - start(eejs.tx[cods.gr[!is.na(cods.gr$eejPrecede)]$eejPrecede]) + 1
    
    cods.gr$eejPrecedeName <- as.character(NA)
    cods.gr[!is.na(cods.gr$eejPrecede)]$eejPrecedeName <-
      eejs.tx$exon_name[cods.gr[!is.na(cods.gr$eejPrecede)]$eejPrecede]
    
    cods.gr[overlapsAny(cods.gr, eejs.tx)]$eejPrecedeDist <- 0
    cods.gr[overlapsAny(cods.gr, eejs.tx)]$eejFollowDist  <- 0 
    
    cods.gr$eejNearestDist <-
      with(mcols(cods.gr), ifelse(
        test = abs(eejFollowDist) < abs(eejPrecedeDist) & !is.na(eejFollowDist), 
        yes  = eejFollowDist,
        no   = eejPrecedeDist ))
    
    cods.gr <- resize(cods.gr, 3, "start")
    
    ### EEJ overlaps
    cods.gr$is_EEJ <-
      overlapsAny(cods.gr, eejs.tx)
    
    cods.exon.hits <-
      findOverlaps(query = resize(cods.gr,1,"start"), subject = ex_on_tx) # find exons for first nucleotide of codon
    
    cods.gr$exon_name <- as.character(NA)
    cods.gr[queryHits(cods.exon.hits)]$exon_name <-
      ex_on_tx$exon_name[subjectHits(cods.exon.hits)]
    
    ### annotate with exon name
    cods.exon.hits <-
      findOverlaps(query = resize(cods.gr,1,"start"), subject = ex_on_tx) # find exons for first nucleotide of codon
    
    cods.gr$exon_name <- as.character(NA)
    cods.gr[queryHits(cods.exon.hits)]$exon_name <-
      ex_on_tx$exon_name[subjectHits(cods.exon.hits)]
    
    ### annotate with CDS exon length (make sure to use the CDS-truncated lengths)
    cods.gr$exon_length <-
      ex_on_tx$exlen[match(cods.gr$exon_name, ex_on_tx$exon_name)]
    
    ### wobble A
    cods.gr$is_wobbleA <-
      grepl("..A", cods.gr$codon)
    
    ### GC3
    cods.gr$is_gc3 <-
      grepl("..[GC]", cods.gr$codon)
    
    return(cods.gr)  
  }


get_codon_matched_sites <-
  function(cods.gr = NULL, m6A_var = "is_m6A"){
    require(doParallel)
    require(GenomicRanges)
    require(hiAnnotator)
    
    ### define codon matched m6A and control sites
    cods.gr.nearest.drach <-
      foreach(cdn = unique(cods.gr[mcols(cods.gr)[[m6A_var]]]$codon), .combine = c, .packages = c("GenomicRanges", "hiAnnotator") ) %dopar% {
        m6A =  cods.gr[  which(   mcols(cods.gr)[[m6A_var]] & cods.gr$codon == cdn) ]
        unm =  cods.gr[  which( ! mcols(cods.gr)[[m6A_var]] & cods.gr$codon == cdn & cods.gr$is_drach) ]
        if(!any(unique(seqnames(m6A)) %in% unique(seqnames(unm)))){
          GRanges(NULL)
        }else{
          #unm[nearest(x = m6A, subject = unm)]
          getNearestFeature(sites.rd =  m6A,  features.rd = unm, colnam = "nearest_unm", feature.colnam = "codonID", parallel = F)
        }
      }
    
  
    cods.gr.nearest.nodrach <-
      foreach(cdn = unique(cods.gr[mcols(cods.gr)[[m6A_var]]]$codon), .combine = c, .packages = c("GenomicRanges", "hiAnnotator") ) %dopar% {
        m6A =  cods.gr[ which(   mcols(cods.gr)[[m6A_var]] & cods.gr$codon == cdn) ]
        unm =  cods.gr[ which( ! mcols(cods.gr)[[m6A_var]] & cods.gr$codon == cdn & ! cods.gr$is_drach) ]
        if(!any(unique(seqnames(m6A)) %in% unique(seqnames(unm)))){
          GRanges(NULL)
        }else{
          #unm[nearest(x = m6A, subject = unm)]
          getNearestFeature(sites.rd =  m6A,  features.rd = unm, colnam = "nearest_unm", feature.colnam = "codonID", parallel = F)
        }
      }
    
    
    is_matched_DRACH <-
      cods.gr$codonID %in% cods.gr.nearest.drach$nearest_unm
    is_matched_noDRACH <-
      cods.gr$codonID %in% cods.gr.nearest.nodrach$nearest_unm
    
    matched_status <-
      ifelse(test = mcols(cods.gr)[[m6A_var]],
             yes  = "m6A",
             no   = ifelse(test = is_matched_DRACH,
                           yes  = "DRACH",
                           no   = ifelse(test = is_matched_noDRACH,
                                         yes  = "no DRACH",
                                         no   = NA)))
    
    matched_status <-
      factor(matched_status, levels=c("no DRACH", "DRACH", "m6A"))
    
    return(matched_status)
    
  }

get_drach_status <-
  function(cods.gr = NULL, m6A_var = "is_m6A"){
    drach_codons <- c("AAA", "GAA", "TAA", "AGA", "GGA", "TGA", "GAC", "AAC", "ACA", "ACC", "ACT")
    stop_codons <- c("TGA", "TAA", "TAG")
    drach_status <-
      ifelse( test = ! cods.gr$codon %in% setdiff(drach_codons, stop_codons),
                                   yes = NA,
                                   no = ifelse(test = mcols(cods.gr)[[m6A_var]],
                                               yes = "m6A",
                                               no = ifelse(test = cods.gr$is_drach,
                                                           yes = "DRACH",
                                                           no = "no DRACH")))
    drach_status = factor(drach_status, levels=c("no DRACH", "DRACH", "m6A"))
    return(drach_status)
  }

get_kmer <- ## get sequence of a small kmer around the site
  function(gr, txs = tx.seqs){
    kmer <-
      rep(NA, length(gr))
    in_txs <-
      as.logical(seqnames(gr) %in% names(txs))
    in_gr <-
      as.logical(names(txs) %in% seqlevels(gr))
    kmer[in_txs] <-
      as.character(mygetSeq(txs[in_gr], gr[in_txs]))
    return(kmer)
  }

get_region <- ## annotate with transcript region (UTR5, CDS, UTR3)
  function(gr, txmodel = txdat$txmodel){
  require(data.table)
    ifelse(
       test = start(gr) < txmodel[match(as.character(seqnames(gr)), tx_name)]$cds.start , # beware of 1-based vs. 0-based
       yes  = "utr5",
       no   = ifelse(
         test = start(gr) >= txmodel[match(as.character(seqnames(gr)), tx_name)]$cds.start & 
                start(gr)  < txmodel[match(as.character(seqnames(gr)), tx_name)]$utr3.start,
         yes  = "cds",
         no   = ifelse(
           test = start(gr) >= txmodel[match(as.character(seqnames(gr)), tx_name)]$utr3.start,
           yes  = "utr3",
           no   = NA
        )
       )
      )
  }


get_region_exex <- ## annotate with transcript region (UTR5, CDS, UTR3)
  function(gr , txmodel = txdat$txmodel){
    require(data.table)
    reg <-
    ifelse(
      test = start(gr) < txmodel[match(as.character(seqnames(gr)), tx_name)]$cds.start , # beware of 1-based vs. 0-based
      yes  = "utr5",
      no   = ifelse(
        test = start(gr) >= txmodel[match(as.character(seqnames(gr)), tx_name)]$utr3.start,
        yes  = "utr3",
        no   = ifelse(
            test = start(gr) >= txmodel[match(as.character(seqnames(gr)), tx_name)]$cds.start & 
                   start(gr)  < txmodel[match(as.character(seqnames(gr)), tx_name)]$utr3.start &
                   start(gr)  < txmodel[match(as.character(seqnames(gr)), tx_name)]$exex_pos_last ,
            yes  = "cds.before",
            no   = ifelse(
              test = start(gr) >= txmodel[match(as.character(seqnames(gr)), tx_name)]$cds.start & 
                     start(gr)  < txmodel[match(as.character(seqnames(gr)), tx_name)]$utr3.start &
                     start(gr) >= txmodel[match(as.character(seqnames(gr)), tx_name)]$exex_pos_last ,
              yes  = "cds.after",
              no   = NA
            )
        )
      )
    )
    reg <- ifelse(is.na(txmodel[match(as.character(seqnames(gr)), tx_name)]$exex_pos_last) & !grepl("utr", reg) , "cds.single", reg )
    reg <- factor(reg, levels = c("utr5", "cds.single", "cds.before", "cds.after", "utr3")) 
    return(reg)
  }



## annotate exon junctions:
## "single" for tx with only one exon-exon junction
## "first"  for tx the first exex
## "internal"  for internal exex (if there are more than two)
## "last"   for the final exex




## ----- functions for generating the windows ------

get_cds_wins <-
  function(txmodel, winsize){
    start <-
      txmodel[!is.na(cds.start), GRanges(seqnames = tx_name, ranges = IRanges(start=cds.start, width=1))] 
    
    stop  <-
      txmodel[!is.na(cds.start), GRanges(seqnames = tx_name, ranges = IRanges(start=cds.start+ cds_len -3, width=1))]
   
  
  cds <- list(
    start = win_filter(start + winsize, txmodel),
    stop  = win_filter(stop  + winsize,  txmodel)
  )
  
  for(i in seq_along(cds)){
    names(cds[[i]])   =   paste0(names(cds)[i], "_", (1:length(cds[[i]])))
  }
  return(cds)
}


get_exex_wins <- 
  function(txmodel = txdat$txmodel, winsize){
    require(data.table)
    # subset to only transcripts with exex
    tx <- txmodel[nexon > 1, ]
    # create GRanges
    gr <- 
      GRanges(seqnames = rep(tx$tx_name, unlist(lapply(tx$exex_anno, length)) ),
              ranges   = unlist(as(sapply(tx$exex_pos, FUN = function(x){IRanges(start=x,end=x)}), "IRangesList")),
              strand ="+"
      ) #+ winsize    
    # fetch annotation
    an <- 
      unlist(tx[,exex_anno]) 
    # group GRanges by annotation
    grl <- list() 
    
    for(i in unique(an)){
      grl[[i]] <- gr[an == i]
      
    }
    
    grl <- lapply(grl, FUN = function(x){win_filter(x + winsize)})
    
    names(grl) <- paste0("exex.",unique(an))
    for(i in seq_along(grl)){ names(grl[[i]])   =   paste0(names(grl)[i], "_", (1:length(grl[[i]])) )}
    return(grl)
  }


get_m6A_wins <-
  function (m6Afile, winsize, txmodel = txdat$txmodel, ex_by_tx = txdat$ex_by_tx,  txs = tx.seqs) {
    ### make m6A windows in transcriptome space
    ###
    bed <-
      read.delim(
        m6Afile,
        sep = "\t",
        header = F,
        stringsAsFactors = F
      )
    colnames(bed) <- c('chr', 'start', 'end', 'id', 'score', 'strand')
    
    bed.gr <- 
      makeGRangesFromDataFrame(
        bed,
        keep.extra.columns =T,
        starts.in.df.are.0based = T,
        seqinfo = Seqinfo(seqnames = as.character(unique(bed$chr)))
      )
    #names(bed.gr) <- mcols(bed.gr)$id
    
    m6A <- mapToTranscripts(bed.gr, ex_by_tx)
    strand(m6A) <- "+" # in tx space, only + strand is allowed. this is not set by mapToTranscripts!
    mcols(m6A) <- NULL # we need to remove mcols in order to merge all windows in one GRanges object
    ### annotate with some metadata
    
    # mcols(m6A)$region <-
    #   get_region(m6A, txmodel)
    # mcols(m6A)$frame <-
    #   get_frame(m6A, txmodel)
    # mcols(m6A)$codon <-
    #   get_codon(m6A, txmodel, txs)
    
    m6A.win <- win_filter(m6A + winsize)
    
    ## get all DRACH in tx with m6A win
    drach <- GRanges(vmatchPattern("DRACH", txs[names(txs) %in% seqlevelsInUse(m6A)], fixed = F)) -2 # center on A
    ## remove known m6A sites
    drach <- drach[!overlapsAny( drach, m6A + winsize )]
    
    # mcols(drach)$region <-
    #   get_region(drach, txmodel)
    # 
    # mcols(drach)$frame <-
    #   get_frame(drach, txmodel)
    # mcols(drach)$codon <-
    #   get_codon(drach, txmodel, txs)
    
    ###
    ###
    drach.win <- win_filter(drach + winsize)
    
    grl <- 
      list(m6A = m6A.win,
           unm = drach.win)
    
    for(i in seq_along(grl)){
      names(grl[[i]])   =   paste0(names(grl)[i], "_", (1:length(grl[[i]])))
    }
    
    return(grl)
    
  }


get_codon_wins <-
  function(txmodel = NULL, winsize = NULL){
    ## make a GRanges of CDS in transcriptome space
    cds.gr <-
      makeGRangesFromDataFrame(
        df = txmodel[!is.na(cds.start),], # the tx_subs makes sure we only use valid transcripts
        seqnames.field = "tx_name",
        start.field = "cds.start",
        end.field = "utr3.start",
        starts.in.df.are.0based = F # the tx model is 1-based
      ) 
    ## resize to exclude first nucleotide of 3'UTR (comes from using 1-based utr3.start as end.field for creating the GRanges)
    cds.gr <- resize(x = cds.gr, width = width(cds.gr)-1, fix="start")
    
    cds.gr <- cds.gr[(width(cds.gr) ) %% 3 == 0] # additional filter to exclude CDSs that are not a multiple of 3 
    cds.gr <- cds.gr [seqnames(cds.gr) %in% names(tx.seqs)] # also, we remove CDSs that are not annotated by sequence (sometimes happens with custom genomes) because this breaks getSeq
    
    ##tile into individual codons  
    cods.gr <-
      unlist(
        tile(
          cds.gr,
          width = 3
        )
      )
    
    cods.gr <- resize(cods.gr, 1, fix = "start") # center on first nt of codon
    names(cods.gr) <-
      paste0("cdn_", 1:length(cods.gr))
    
    grl <- list(all_codons = win_filter(cods.gr + winsize))
    return(grl)
    
    
  }


get_kmer_wins <-
  
  function(kmer, winsize, txmodel = txdat$txmodel, ex_by_tx = txdat$ex_by_tx,  txs = tx.seqs) {
    
    kmer.gr  <- GRanges(vmatchPattern(kmer, txs, fixed = F)) -2 # center on A
    kmer.win <- win_filter(kmer.gr + winsize)
     
    grl <-
      list(kmer.win)
    names(grl) = kmer

    for(i in seq_along(grl)){
      names(grl[[i]])   =   paste0(names(grl)[i], "_", (1:length(grl[[i]])))
    }

    return(grl)
    
  }

