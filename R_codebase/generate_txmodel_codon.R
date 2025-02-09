###
### ## the coordinates of the features in the txmodel table are 1-based!!!
###

### ----- helper functions -----

## get exon-exon junction positions in tx
get_exex_pos <- function(txmodel, ex_by_tx){
  exex_pos <- cumsum( width(ex_by_tx[match(txmodel$tx_name,names(ex_by_tx))]) ) +1 # !!! make 1-based !!! to be consistent with the rest of txmodel
  exex_pos <- lapply(exex_pos, FUN = function(x) unlist(x))           # make list of vectors parallel to txmodel 
  exex_pos[lapply(exex_pos, FUN = function(x) length(x)) <=1] <- NA    # single exon genes do not have exex junction
  exex_pos <- lapply(exex_pos, FUN = function(x) x[1:(length(x)-1)] )  # the end of the last exon is not an exex junction
  return(exex_pos)  
}

### function for annotating exon-exon junction into categories:
## single:  only one exon junction in tx (i.e. two exon tx)
## first:   first exon junction in transcript
## last:    last exon junction in transcript
## this function returns a list of vectors that is parallel to the exex_pos list/column in the txmodel
## use lapply to apply it in the exex_pos column
exex_anno_vec <- 
  function(x){
    if( any(is.na(x) | length(x)==0 ) ){NA}
    else if (length(x) == 1){"single"}
    else{ifelse(test = x == min(x),
                yes  = "first",
                no   = ifelse(test = x == max(x),
                              yes  = "last",
                              no   = "internal"
                )
    )
    }
  }



get_gene_flanks <- function(grl, fastafile,  flank_up, flank_down){
  require(Biostrings)
  chrom_len <- seqlengths(readDNAStringSet(fastafile))
  extab <- data.table(tx_name = names(unlist(grl)),
                      exon_rank = grl@unlistData@elementMetadata@listData$exon_rank)
  extab <- extab[ , .(ex_last = max(exon_rank)), by = tx_name ]
  
  singleex <- unlist(grl)
  
  mcols(singleex)$firstex = mcols(singleex)$exon_rank == 1
  mcols(singleex)$lastex  = mcols(singleex)$exon_rank == extab[match(names(singleex), extab$tx_name),]$ex_last
  
  singleex[mcols(singleex)$firstex, ] <-
    resize(singleex[mcols(singleex)$firstex, ],
           width = width(singleex[mcols(singleex)$firstex, ]) + flank_up,
           fix = "end")
  singleex[mcols(singleex)$lastex,] <-
    resize(singleex[mcols(singleex)$lastex,],
           width = width(singleex[mcols(singleex)$lastex, ])  + flank_down,
           fix = "start")
  
  ## generate filter for out of bounds GRanges
  mcols(singleex)$oob <- FALSE
  mcols(singleex[start(singleex) < 0])$oob <- TRUE
  mcols(singleex[end(singleex) > 
                   chrom_len[as.integer(match(seqnames(singleex), names(chrom_len)))]])$oob <- TRUE
  
  #singleex[!mcols(singleex)$oob]
  
  grl <- relist(singleex, skeleton = grl)
  
  return(grl)
  
}


### ------ main function to call ------

get_txdat_from_gtf <- function(txfile, flank_up = 0, flank_down = 0, fastafile = NULL){
  require("GenomicFeatures")
  require("data.table")
  
  if( flank_up != 0 | flank_down != 0  ){
    if(is.null(fastafile)){
      stop("Flank extension requires a valid genome fasta file.")
    }else{if(!file.exists(fastafile)){
        stop("Flank extension requires a valid genome fasta file.")
        }
      }
    }
  

   # read in gtf
   txdb <- makeTxDbFromGFF(txfile)
   # generate list of exons grouped by tx
   ex_by_tx <-
    exonsBy(txdb, by = "tx", use.names = TRUE) ### this gives positions of exex per tx
   
   ### here, we modify the gene flanks ### This is work in progress and buggy!!!
   if(flank_up != 0 | flank_down != 0){
     ex_by_tx <- get_gene_flanks(grl = ex_by_tx, flank_up = flank_up, flank_down = flank_down, fastafile = fastafile)
     
     oob <-     data.table(tx_name = do.call(rbind, strsplit(names(unlist(ex_by_tx)), split = "\\.") )[,1] ,
                           oob = mcols(unlist(ex_by_tx))$oob)
     has.oob <- oob[, oob, by = tx_name]
     
     ex_by_tx <- ex_by_tx[!has.oob[match(names(ex_by_tx), has.oob$tx_name),]$oob]
     
     }
   
   
  # generate data table describing tx features
  txmodel <-
    data.table(transcriptLengths(
      txdb,
      with.utr5_len = T,
      with.cds_len = T,
      with.utr3_len = T
    )) ### this gives positions of cds start and stop per tx
  
  txmodel <- txmodel[txmodel$tx_name %in% names(ex_by_tx)]
  
  if(flank_up != 0 | flank_down != 0){
  ## the flank feature is needed to extend UTRs, e.g. iun yeast
  txmodel$utr5_len = txmodel$utr5_len + flank_up
  txmodel$utr3_len = txmodel$utr3_len + flank_down
  txmodel$tx_len   = txmodel$tx_len   + flank_up + flank_down
}else{
  ## the flank feature is needed to extend UTRs, e.g. iun yeast
  txmodel$utr5_len = txmodel$utr5_len #+ flank_up
  txmodel$utr3_len = txmodel$utr3_len #+ flank_down
  txmodel$tx_len   = txmodel$tx_len   #+ flank_up + flank_down
}
  ## for transcript regions of length 0, set NA
  txmodel[utr5_len == 0, utr5_len := NA]
  txmodel[cds_len  == 0, cds_len  := NA]
  txmodel[utr3_len == 0, utr3_len := NA]
  
  ## determine region boundaries
  txmodel[, cds.start := utr5_len +1 ] # first nucleotide of start codon !!! 1-based !!!
  #txmodel[, cds.stop  := utr5_len +1 + cds_len  -3] # first nucleotide of stop codon
  txmodel[, utr3.start:= utr5_len +1 + cds_len ] # first nucleotide of utr3
  txmodel[is.na(utr3_len), utr3.start := NA] # we need to remove utr3.start from transcripts without utr3
  
  ## determine exon-exon junction positions
  txmodel$exex_pos <-
    get_exex_pos(txmodel, ex_by_tx)
  
  txmodel$exex_anno <-
    lapply(txmodel$exex_pos, exex_anno_vec)

  ## determine positions of first, penultimate and last exex junctions
  txmodel[nexon  > 2, exex_pos_first  := unlist(lapply(exex_pos, min))]
  txmodel[nexon  > 2, exex_pos_last   := unlist(lapply(exex_pos, max))]
  txmodel[nexon == 2, exex_pos_last := unlist(lapply(exex_pos, max))] ## as we are interested in the distinction between final exon and prior exons, we include the final exon of two exon tx
  
  ## determine length of cds before and after last exex # todo: double-check for 1-based calculations
  txmodel[exex_pos_last > cds.start & exex_pos_last < utr3.start, # we need to filter exex in utrs
          cds.before_len := exex_pos_last - cds.start]
  txmodel[exex_pos_last > cds.start & exex_pos_last < utr3.start, # we need to filter exex in utrs
          cds.after_len  := utr3.start    - exex_pos_last]
  
  ## Determine exon junction density in CDS
 
  # require(GenomicRanges)
  # require(GenomicFeatures)
  # ## relative position of last exex in cds
  # exex_last_pos_rel.cds    <-  (txmodel$exex_pos_last - txmodel$cds.start)  / txmodel$cds_len
  # exex_last_pos_rel.cds[exex_last_pos_rel.cds < 0 | exex_last_pos_rel.cds > 1] <- NA
  # txmodel$exex_last_pos_rel.cds <- exex_last_pos_rel.cds
  # 
  # ## position of penultimate exex junction
  # exex_pos_pen <- lapply(txmodel$exex_pos, FUN = function(x) x[order(x, decreasing = TRUE) == 2] )
  # exex_pos_pen <- as.numeric(exex_pos_pen)
  # txmodel$exex_pos_pen <- exex_pos_pen
  # 
  # ## exon lengths
  # exon_len <- width(ex_by_tx[match(txmodel$tx_name,names(ex_by_tx))])
  # exon_len <- lapply(exon_len, FUN = function(x) unlist(x))
  # txmodel$exon_len <- exon_len
  # 
  # ## length of longest exon
  # txmodel$exon_len_max <- unlist(lapply(txmodel$exon_len, max))
  
  # ## intron lengths
  # int_by_tx <- intronsByTranscript(txdb)
  # intron_len <- width(int_by_tx[match(txmodel$tx_name,names(ex_by_tx))])
  # intron_len <- lapply(intron_len, FUN = function(x) unlist(x))
  # txmodel$intron_len <- intron_len
  
  ## determine other features like exex in cds and utrs
  
  exex_in_utr5 <- function(onerow)  { lapply(onerow$exex_pos, FUN = function(x) x < onerow$cds.start )} 
  exex_utr5    <-  apply(txmodel,1,exex_in_utr5 )
  exex_utr5    <-  lapply(exex_utr5, FUN = function(x) sum(unlist(x))  ) 
  txmodel$exex_utr5 <- unlist(exex_utr5)
  
  
  exex_in_cds <- function(onerow)  { lapply(onerow$exex_pos, FUN = function(x) x > onerow$cds.start &  x < onerow$utr3.start )} 
  exex_cds    <-  apply(txmodel,1,exex_in_cds )
  exex_cds    <-  lapply(exex_cds, FUN = function(x) sum(unlist(x))  ) 
  txmodel$exex_cds   <- unlist(exex_cds)
  txmodel$exex_cds_dens <-  txmodel$exex_cds / txmodel$cds_len  
  
 
  exex_in_utr3 <- function(onerow)  { lapply(onerow$exex_pos, FUN = function(x) x > onerow$utr3.start )} 
  exex_utr3    <-  apply(txmodel,1,exex_in_utr3 )
  exex_utr3    <-  lapply(exex_utr3, FUN = function(x) sum(unlist(x))  ) 
  txmodel$exex_utr3   <- unlist(exex_utr3)
  
  ###
  ### novel method to determine exex cds density / exon length 
  ###
  get_cds_exon_len <- function(onerow){
    
    if(onerow$nexon == 1){
      onerow$cds_len
    }else if(is.na(onerow$cds.start)|is.na(onerow$utr3.start)){NA}else{
      
      mod.start <-
        c(onerow$cds.start,

          unlist(onerow$exex_pos)[ unlist(onerow$exex_pos) >  onerow$cds.start &
                                     unlist(onerow$exex_pos) <  onerow$utr3.start]

        )
      mod.end <- c(unlist(onerow$exex_pos)[ unlist(onerow$exex_pos) >  onerow$cds.start &
                                              unlist(onerow$exex_pos) <  onerow$utr3.start ],
                   onerow$utr3.start ) -1
      width(
        IRanges(
          start=mod.start,
          end=mod.end)
      )

    } 
    
  }
  
  txmodel$cds_exon_len <-
    apply(txmodel,1,get_cds_exon_len )
  
  ### create idtab
  gtf <-
      unique(as.data.table(mcols(rtracklayer::import(txFile))))
    gtf[, ensembl_transcript_id := gsub("\\..*$", "", transcript_id)]
    gtf[, ensembl_gene_id := gsub("\\..*$", "", gene_id)]
    gtf[, external_gene_name := gsub("\\..*$", "", gene_name)]
    
    idtab <-
      gtf[, c("ensembl_gene_id",
              "external_gene_name",
              "ensembl_transcript_id", 
              "exon_id")]
    idtab <-
      idtab[!is.na(ensembl_transcript_id),]
    idtab <-
      unique(idtab)
  
  ## remove duplicate external_gene_names
  if(length(unique(unique(idtab[,c("ensembl_gene_id", "external_gene_name") ])$ensembl_gene_id)) !=
     length(unique(unique(idtab[,c("ensembl_gene_id", "external_gene_name") ])$external_gene_name))){
    
    rbind(unique(unique(idtab[,c("ensembl_gene_id", "external_gene_name") ]))[!isUnique(ensembl_gene_id),][order(ensembl_gene_id),][, error := "Non-unique ID"],
          unique(unique(idtab[,c("ensembl_gene_id", "external_gene_name") ]))[!isUnique(external_gene_name),][order(external_gene_name),][, error := "Non-unique Name"])
    
    idtab <-
      idtab[ ! ensembl_gene_id %in% unique(idtab[, c("ensembl_gene_id", "external_gene_name"), with=F])[!isUnique(external_gene_name), ][, .SD[1], by=c("external_gene_name")]$ensembl_gene_id, ]
    
    warning("Non-unique ID-to-Name relations detected and removed!\n")
    
  }
  
    get_intronic_parts <-
      ## extract exonic and intronic regions per gene to enable separate read counting
      ## adapted from: https://gist.github.com/LTLA/6fde8edf558a0d32314043ab59e253cf#file-intron_exon_counter-r
      function(txdb = NULL){
        introns.gr <-
          intronicParts(txdb, linked.to.single.gene.only=TRUE)
        introns.dt <-
          data.table(
            GeneID = introns.gr$gene_id,
            Chr = as.character(seqnames(introns.gr)),
            Start = as.integer(start(introns.gr)),
            End = as.integer(end(introns.gr)),
            Strand = as.character(strand(introns.gr)),
            stringsAsFactors = FALSE
          )
        
        exons.gr <-
          exonicParts(txdb, linked.to.single.gene.only = TRUE)
        exons.dt <-
          data.table(
            GeneID = exons.gr$gene_id,
            Chr = as.character(seqnames(exons.gr)),
            Start = as.integer(start(exons.gr)),
            End = as.integer(end(exons.gr)),
            Strand = as.character(strand(exons.gr)),
            stringsAsFactors = FALSE
          )
        return(list(exons.gr = exons.gr, exons.dt = exons.dt, introns.gr = introns.gr, introns.dt = introns.dt))
      }
    intronic_parts <-
      get_intronic_parts(txdb)
 
     
    
    
  ### combine into txdat object
  txdat <- list( txmodel =  txmodel, ex_by_tx = ex_by_tx, idtab = idtab, intronic_parts = intronic_parts)
  
  return(txdat)
 
}



