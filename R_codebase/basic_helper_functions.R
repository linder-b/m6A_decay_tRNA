#' Make a BSgenome from a fasta file.
#' If required,  a2bit file and the R package will be built and installed
#' 
#' @param fastaFile A path to a fasta file.
#' @param organism Name of organism; used for pkg name.
#' @param genome Name of genome build; used for pkg name.
#' @returns A BSgenome object.
#' @examples
#' create_BSgenome("gencode.v35.transcripts.fa.gz")
create_BSgenome <-
  function(fastaFile = NULL, organism = NULL, genome = NULL){
    require(rtracklayer)
    require(BSgenomeForge)
    
    BSDir = file.path(dirname(fastaFile), "BSgenome")
    if(!dir.exists(BSDir))dir.create(BSDir)
    
    twobitFile = file.path(BSDir, paste0(basename(fastaFile), ".2bit"))
    
    if(!file.exists(twobitFile)){
      cat(paste("No 2bit file found at", twobitFile, ", generating it..."))
      rtracklayer::export.2bit(object = replaceAmbiguities(readDNAStringSet(filepath = fastaFile)), con = twobitFile )
      cat(paste("...done!\n"))
    }else{
      cat(paste("Loading 2bit file found at", twobitFile, "\n"))
    }
    
    
    provider = "Custom"
    organism = BSgenomeForge:::format_organism(organism)
    abbr_organism = BSgenomeForge:::abbreviate_organism_name(organism)
    pkgname = BSgenomeForge:::make_pkgname(abbr_organism, provider, genome)
    
    
    if(!dir.exists(file.path(BSDir, pkgname) )){
      
      cat(paste("Package not found at", file.path(BSDir, pkgname), ", generating it..."))
      
      forgeBSgenomeDataPkgFromTwobitFile(
        filepath = twobitFile
        , organism = organism
        , provider = provider, 
        , pkg_version = "0.0.1"
        , genome = genome
        , pkg_maintainer = "Bastian Linder <linder@umlaut.bio>"
        , destdir = BSDir 
      )
      
      cat(paste("...done!\n"))
      
    }else{
      cat(paste("Package found at", file.path(BSDir, pkgname), "\n"))
    }
    
    # "/home/linder/genomes/homo_sapiens_GRCh38_v35/BSgenome/BSgenome.Hsapiens.Custom.hg38"
    
    if(!pkgname %in% installed.packages()[, "Package"] ){
      install.packages(pkgs = file.path(BSDir, pkgname), repos = NULL, type = "source")
    }
    
    library(pkgname, character.only=TRUE)
    
    bsg = get(pkgname)
    return(bsg)
    
  }


get_frame <- ## get the reading frame of the first nucleotide of a GRanges in transcriptome space
  function(gr, txmodel = NULL){
    pos <- start(gr)   # beware of 1-based vs. 0-based
    ## need to filter sites outside of CDS, as the modulo approach will otherwise assign a frame
    pos[pos <  txmodel[as.integer(match(seqnames(gr), txmodel$tx_name)),]$cds.start]   <- NA
    pos[pos >= txmodel[as.integer(match(seqnames(gr), txmodel$tx_name)),]$utr3.start ] <- NA
    ## assign frame using modulo of position in CDS
    frame <-   (pos - txmodel[as.integer(match(seqnames(gr), txmodel$tx_name)),]$cds.start) %% 3 # beware of 1-based vs. 0-based
    return(factor(frame))
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

relpos.cds.gr <- function(gr, txmodel = txdat$txmodel){
  require("GenomicAlignments")
  require("GenomicFeatures")
  require("data.table")
  
  dt <-
    data.table(pos  = start(gr),
               txid = as.character(seqnames(gr)))
  
  
  get_relpos <- function(txid, pos, txmodel){
    idx = match(txid, txmodel$tx_name)
    pos.utr5 =  pos / txmodel$utr5_len[idx]
    pos.utr5[pos.utr5 < 0 | pos.utr5 > 1] <- NA
    pos.cds  = (pos - txmodel$cds.start[idx]) / txmodel$cds_len[idx]
    pos.cds[pos.cds < 0 | pos.cds > 1] <- NA
    pos.utr3 = (pos - txmodel$utr3.start[idx]) / txmodel$utr3_len[idx]
    pos.utr3[pos.utr3 < 0 | pos.utr3 > 1] <- NA
    return(list(pos.utr5, pos.cds, pos.utr3))
  }
  
  
  dt[,  c("pos.utr5", "pos.cds", "pos.utr3")  := get_relpos(txid, pos, txmodel)   ][,]
  
  dt[, c("pos.utr5", "pos.utr3") := list(pos.utr5 * (median(txmodel$utr5_len, na.rm=T)/median(txmodel$cds_len, na.rm=T)),
                                         pos.utr3 * (median(txmodel$utr3_len, na.rm=T)/median(txmodel$cds_len, na.rm=T)) )]
  
  
  dt[, c("pos.utr5", "pos.cds", "pos.utr3")  :=  list(pos.utr5 + (1-max(pos.utr5, na.rm=T)),
                                                      pos.cds  + 1,
                                                      pos.utr3 + 2)   ]
  
  dt[ !is.na(pos.utr5), relpos := pos.utr5][ !is.na(pos.cds), relpos := pos.cds][ !is.na(pos.utr3), relpos := pos.utr3]
  
  return(dt$relpos)
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
      #get_kmer(gr = cods.gr, txs = tx.seqs)
      as.character(getSeq(tx.seqs, cods.gr))
    
    cods.gr$kmer <-
      #get_kmer(cods.gr+2, tx.seqs)
      as.character(getSeq(tx.seqs, cods.gr+3))
    
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

## a function that extracts proper gene symbols and removes genes not found in annotation
translate_symbol <-
  function(gene = NULL, dict = NULL, fill = FALSE) {
    dict = as.character(dict)
    dict = unique(dict)
    symbol = dict[match(toupper(gene), toupper(dict))]
    names(symbol) = NULL
    ignored = gene[is.na(symbol)]
    if(fill != TRUE)symbol = symbol[!is.na(symbol)]
    if (length(ignored > 0))
      warning(paste(
        "The following genes were not found and ignored:",
        paste(ignored, collapse = ", "),
        "\n"
      ))
    return(symbol)
  }

### pre-defined genesets
get_geneset <-
  function(gsList = NULL, dict = NULL) {
    ## if no custom list is given we use the default one
    if(is.null(gsList)){
      ## maually curated genesets of interest
      require(ProliferativeIndex)
      gsList = list(
        m6A     = c("METTL3", "METTL14", "WTAP", "RBM15", "ZC3H13", "HAKAI", "CBLL1", "VIRMA", "KIAA1429"),
        m6A_neg = c("ZFP217", "ZNF217", "ALKBH5"),
        yth     = c("YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2"),
        emt     = c("Zeb1", "Zeb2", "Snai1", "Snai2", "Tnc", "Fn1", "Ncam1", "Cdh1", "Cdh2", "Vim", "TGFB", "TWIST1"),
        emt_E   = c("CDH1", "DSP", "EPCAM", "JUP", "KRT18", "KRT19", "KRT7", "KRT8"), ## from Ebright 2020 Science, Fig.S11
        emt_M   = c("SERPINE1", "VIM", "SNAI1", "ZEB1", "ZEB2", "CDH2", "FOXC2",  "SNAI2", "TWIST1", "GSC", "FN1", "ITBG6", "MMP2", "MMP3", "MMP9", "SOX10"),
        igf     = c(paste0("IGF2BP", 1:3)),
        mcm5    = c("IKBKAP", paste0("ELP", 1:6), "ALKBH8"),
        s2U     = c(paste0("CTU", 1:2), "UBA4", "MOCS3", "URM1"), 
        t6A     = c("OSGEP", "LAGE3", "TPRKB", "TPRK", "TP53RK", "TRP53RK", "TRP53RKA", "TRP53RKB", "C14orf142", "GON7"),
        tRNAlig = c("RTCB", "ZBTB8OS", "DDX1", "C2orf49", "FAM98B"),
        ngd     = c("ASCC2", "ASCC3", "RNF25", "TRIP4"),  ## as determined by depmap clustering 
        pluri   = c( "Esrrb", "Dppa3", "Zfp42", "Jarid2", "Klf2", "Nr0b1", "Fgf4", "Nanog"),
        oskm    = c("Pou5f1", "Oct4","Sox2", "Klf4", "Myc"),
        immune  = c("CD274", "CX3CL1", "PDCD1", "CD80", "CTLA4"), ## from Herbst, et al. (2014). Predictive correlates of response to the anti-PD-L1 antibody MPDL3280A in cancer patients. Nature 515, 563–567.
        proliferation = unique(ProliferativeIndex:::metaPCNA2), ## non-unique genes present: ProliferativeIndex:::metaPCNA2[!isUnique(ProliferativeIndex:::metaPCNA2)]
        custom  = c("TP53", "BRAF", "KRAS", "NRAS", "HRAS", "SETD2", "VHL")
      )
    }
    ## if dictionary is given, we adapt gene names
    if(!is.null(dict)){
      gsList <-
        lapply(gsList, function(x)translate_symbol(gene = x, dict = dict))
    }
    print(unlist(lapply(gsList, length)))
    return(gsList)
  }

get_geneset_name <-
  function(gene, geneset){setkeyv(rbindlist(lapply(names(geneset), function(geneset_name){data.table(gene = geneset[[geneset_name]], geneset_name = geneset_name) } )), "gene")[gene, geneset_name]}


## a function to transpose a data table, using colnames as new rownames and the first column values as new colnames
tdt <-
  function(dt = NULL){
    as.data.table(t(as.matrix(dt[, -1], rownames = dt[[1]])), keep.rownames = TRUE)
  }

scale1 <-
  function(x){as.numeric(scale(x,center=min(x, na.rm = TRUE),scale=diff(range(x, na.rm = TRUE))))}

## a function to convert FPKM to TPM; from: https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/307#issuecomment-478092373
## referring to: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
# FPKMtoTPM <- function(x) {
#   return(exp(log(x) - log(sum(x)) + log(1e6)))
# }
## functions directly from: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate), na.rm = TRUE))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}
