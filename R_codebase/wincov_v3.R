# source("generate_windows_w_codons.R")

ribocov_from_bamfile <-
  function(bamName = NULL, offsets = NULL, rwinsize = 60, threads = 1){
    #stopifnot(exists("ex_by_tx", "txmodel"))
    ## define input
    bamfile = bamfiles[[bamName]]
    if(!file.exists(bamfile))stop(paste("bamflile", bamfile, "not found.\n"))
    if(is.null(offsets))stop(paste("No offsets defined.\n"))
    offsets[, rwin := rwinsize]
    chksm  = digest::digest(offsets[run == bamName, c("run", "use", "offset", "rwin")]) ## this checks over all paramaters that determine coverage
    outFile = paste0(tmpDir, "/", gsub("\\.bam$", paste0(".ribocov.raw", ".offset", chksm,".RData"),basename(bamfile)))
    ofsFile = paste0(tmpDir, "/", gsub("\\.bam$", paste0(".offsets", chksm,".RData"),basename(bamfile)))
    message(paste("Start processing", basename(bamfile), "at", Sys.time(), "\n"))
    if(!file.exists(outFile)){
      message(paste("Generating", basename(outFile), "at", Sys.time(), "\n"))
      ## which transcripts (only those with windows)
      bamwhich = unlist(GenomicRanges::reduce(ex_by_tx))
      seqlevels(bamwhich, pruning.mode="coarse") <- seqlevels(BamFile(bamfile)) ## adapt seqlevels to bam file
      #names(which) <- NULL
      #return(bamwhich)
      ## scanBam parameters
      bamparam.wins <-
        ScanBamParam(
          what = c("qwidth", "cigar"),
          ## qwidth is essential for FP size calculations
          flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE),
          which = bamwhich 
        )
      
      ## determine if paired  
      is_paired <-
        suppressMessages(testPairedEndBam(BamFile(bamfile)))
      ## get GRanges
      if (is_paired) {
        gr <-
          GRanges(readGAlignmentPairs(
            file = bamfile,
            use.names = T,
            param = bamparam.wins
          ))
      } else{
        gr <-
          GRanges(readGAlignments(
            file = bamfile,
            use.names = T,
            param = bamparam.wins
          ))
      }
      ## update seqlevels
      seqlevels(gr, pruning.mode = "coarse") <-
        seqlevels(ex_by_tx)
      ### map read starts to transcripts
      tx <-
        suppressWarnings(mapToTranscripts(resize(gr, 1, "start"), ex_by_tx)) ## lift read starts to tx space
      
      ## update metadata
      mcols(tx) <-
        cbind(mcols(tx), mcols(gr[tx$xHits]))
      
      ## update strand
      strand(tx) <-
        "+"  ## set strand to + as mapToTranscripts retains genomic strand
      
      ## update seqlengths to get full coverage
      seqinfo(tx) <-
        Seqinfo(seqnames =  seqlevels(tx),
                seqlengths = txmodel$tx_len[match(seqlevels(tx), txmodel$tx_name)])
      
      ### filter fragment sizes (based on above offset metrics)
      tx <-
        tx[tx$qwidth %in% offsets[(use), qwidth]]
      ### Shift to A-sites (based on above offset metrics)
      tx <-
        unlist(as(foreach(i = split(
          x = tx,
          f = as.factor(tx$qwidth),
          drop = T
        )) %dopar% {
          suppressWarnings(shift(i, offsets$offset[match(i[1]$qwidth, offsets$qwidth)]))
        }, "GRangesList"), use.names = F)
      
      ### calculate coverage
      ribocov.raw <-
        coverage(tx)
      ### save coverage file
      saveRDS(object = ribocov.raw, file = outFile)
      ### save offset file
      saveRDS(object = offsets, file = ofsFile)
    }else{
      message(paste("Loading", basename(outFile), "at", Sys.time(), "\n"))
      ribocov.raw = readRDS(outFile)}
    return(ribocov.raw)
  }

covAsite <-
  function(bamFile = NULL, covFile = NULL, offLst = NULL, winSize = 60) {
    require("GenomicRanges")
    if (is.null(offLst))
      stop(paste("No offset list defined.\n"))
    if (!file.exists(bamFile))
      stop(paste("BAM file", bamFile, "not found.\n"))
    if (!file.exists(covFile)) {
      ## setup bam parameters
      bamWhich = unlist(GenomicRanges::reduce(ex_by_tx))
      seqlevels(bamWhich, pruning.mode = "coarse") <-
        seqlevels(BamFile(bamFile)) ## adapt seqlevels to bam file
      bamparam.wins <-
        ScanBamParam(
          what = c("qwidth", "cigar"),
          ## qwidth is essential for FP size calculations
          flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE)
        )
      ## determine if paired
      is_paired <-
        suppressMessages(testPairedEndBam(BamFile(bamFile)))
      ## get GRanges
      if (is_paired) {
        gr <-
          GRanges(readGAlignmentPairs(
            file = bamFile,
            use.names = T,
            param = bamparam.wins
          ))
      }else{
        gr <-
          GRanges(readGAlignments(
            file = bamFile,
            use.names = T,
            param = bamparam.wins
          ))
      }
      ## filter sizes
      gr = gr[gr$qwidth %in% as.numeric(names(offLst))]
      ## update seqlevels
      seqlevels(gr, pruning.mode = "coarse") <-
        seqlevels(ex_by_tx)
      ## map read starts to transcripts
      tx = suppressWarnings(mapToTranscripts(resize(gr, 1, "start"), ex_by_tx))
      print(paste(length(gr), "genome-mapped reads resulted in", length(tx), "transcriptome-mapped reads.\n"))
      ## update metadata (qwidth & cigar)
      mcols(tx) = cbind(mcols(tx), mcols(gr[tx$xHits]))
      ## update strand
      strand(tx) = "+"  ## set strand to + as mapToTranscripts retains genomic strand
      ## determine relative frame (also for pos outside CDS)
      tx$frame = (start(tx) - txmodel[as.integer(match(seqnames(tx), txmodel$tx_name)), ]$cds.start) %% 3
      ## update seqlengths from txmodel to get full coverage with cov()
      seqinfo(tx) = Seqinfo(seqnames =  seqlevels(tx), seqlengths = txmodel[match(seqlevels(tx), tx_name), tx_len])
      ## Shift to A-sites (based on above offset metrics)
      txl = split(x = tx,
                  f = as.factor(paste0(tx$qwidth, "_", tx$frame)),
                  drop = T) ## split by qwidth and frame
      txl = foreach(tx = txl) %dopar% {
        offset =  offLst[[as.character(tx[1]$qwidth)]][[as.character(tx[1]$frame)]] ## shift be offset
        suppressWarnings(GenomicRanges::shift(tx, offset))
      }
      tx = unlist(as(txl, "GRangesList"), use.names = F)
      ## transcriptome coverage
      cov = coverage(tx)
      ## save coverage
      saveRDS(object = cov, file = covFile)
    } else{
      ## load coverage
      cov = readRDS(file = covFile)
    }
    return(cov)
  }

fivePcov_from_bamfile <-
  function(bamName = NULL, rwinsize = 60, threads = 1){
    #stopifnot(exists("ex_by_tx", "txmodel"))
    ## define input
    bamfile = bamfiles[[bamName]]
    if(!file.exists(bamfile))stop(paste("bamflile", bamfile, "not found.\n"))
    chksm  = digest::digest(offsets[run == bamName, c("run", "use", "rwin")]) ## this checks over all paramaters that determine coverage
    outFile = paste0(tmpDir, "/", gsub("\\.bam$", paste0(".fivePcov.raw", ".offset", chksm,".RData"),basename(bamfile)))
    message(paste("Start processing", basename(bamfile), "at", Sys.time(), "\n"))
    if(!file.exists(outFile)){
      message(paste("Generating", basename(outFile), "at", Sys.time(), "\n"))
      ## which transcripts (only those with windows)
      bamwhich = unlist(GenomicRanges::reduce(ex_by_tx))
      seqlevels(bamwhich, pruning.mode="coarse") <- seqlevels(BamFile(bamfile)) ## adapt seqlevels to bam file
      #names(which) <- NULL
      #return(bamwhich)
      ## scanBam parameters
      bamparam.wins <-
        ScanBamParam(
          what = c("qwidth", "cigar"),
          ## qwidth is essential for FP size calculations
          flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE),
          which = bamwhich 
        )
      
      ## determine if paired  
      is_paired <-
        suppressMessages(testPairedEndBam(BamFile(bamfile)))
      ## get GRanges
      if (is_paired) {
        gr <-
          GRanges(readGAlignmentPairs(
            file = bamfile,
            use.names = T,
            param = bamparam.wins
          ))
      } else{
        gr <-
          GRanges(readGAlignments(
            file = bamfile,
            use.names = T,
            param = bamparam.wins
          ))
      }
      
      ### map read starts to transcripts
      tx <-
        suppressWarnings(mapToTranscripts(resize(gr, 1, "start"), ex_by_tx)) ## lift read starts to tx space
      
      ## update metadata
      mcols(tx) <-
        cbind(mcols(tx), mcols(gr[tx$xHits]))
      
      ## update strand
      strand(tx) <-
        "+"  ## set strand to + as mapToTranscripts retains genomic strand
      
      ## update seqlengths to get full coverage
      seqinfo(tx) <-
        Seqinfo(seqnames =  seqlevels(tx),
                seqlengths = txmodel$tx_len[match(seqlevels(tx), txmodel$tx_name)])
      
      ### filter fragment sizes (based on above offset metrics)
      tx <-
        tx[tx$qwidth %in% offsets[(use), qwidth]]
      ### calculate coverage
      cov.raw <-
        coverage(tx)
      ### save coverage file
      saveRDS(object = cov.raw, file = outFile)
      }else{
      message(paste("Loading", basename(outFile), "at", Sys.time(), "\n"))
      cov.raw = readRDS(outFile)}
    return(cov.raw)
  }

## a wrapper to generate normalized coverage
ribonorm_from_ribocov <-
  function(bamName = NULL, offsets = NULL, rwinsize = NULL, threads = 1){
    offsets[, rwin := rwinsize]
    chksm  = digest::digest(offsets[run == bamName, c("run", "use", "offset", "rwin")]) ## this checks over all paramaters that determine coverage
    inFile  = paste0(tmpDir, "/", gsub("\\.bam$", paste0(".ribocov.raw",  ".offset", chksm,".RData"),basename(bamfiles[bamName])))
    outFile = paste0(tmpDir, "/", gsub("\\.bam$", paste0(".ribocov.norm", ".offset", chksm,".RData"),basename(bamfiles[bamName])))
    if(!file.exists(inFile))stop(paste("Input file", inFile, "not found!\n"))
    if(!file.exists(outFile)){
      covnorm = scalecov(cov.raw = readRDS(inFile), rwinsize = rwinsize, threads = threads)
      saveRDS(object = covnorm, file = outFile)
      }else{
      covnorm = readRDS(file = outFile)
      }
    return(covnorm)
    }


covNorm <-
  function(covFile = NULL, covNormFile = NULL, rwinsize = 60, threads = 1){
    if(!file.exists(covFile))stop(paste("Input file", covFile, "not found!\n"))
    if(!file.exists(covNormFile)){
      covNorm = scalecov(cov.raw = readRDS(covFile), rwinsize = rwinsize, threads = threads)
      saveRDS(object = covNorm, file = covNormFile)
    }else{
      covNorm = readRDS(file = covNormFile)
    }
    return(covNorm)
  }

### normalize coverage
#### scaling function for zero variance instead of NA from https://stackoverflow.com/a/15364319
myscale <-
  function(x) {
    (x - mean(x)) / sd(x) ^ as.logical(sd(x))
  }

#### function to apply on coverage object
scalecov <-
  function(cov.raw = NULL, rwinsize = NULL, threads = 1) {
    as(
      mclapply(
        cov.raw,
        FUN = function(tx) {
          rollapply(
            data = as.numeric(tx),
            width = (2 * rwinsize) + 1,
            fill = NA,
            by = 1,
            FUN = function(x) {
              myscale(x)[rwinsize + 1]
            }
          )
        },
        mc.cores = threads,
        mc.cleanup = TRUE
      ),
      "SimpleRleList"
    )
  }


ribodat_from_bamfile <- function(bamfile, txmodel) {
  
  ## random sample of start and stop codons
  set.seed(42)
  wins <-
    unlist(GRangesList(get_cds_wins(txmodel = txmodel[sample(nrow(txmodel), 10000),], winsize = rwinsize)), use.names = FALSE)
  #return(wins)
  ## which transcripts (only those with windows)
  bamwhich = unlist(GenomicRanges::reduce(ex_by_tx[names(ex_by_tx) %in% seqlevelsInUse(wins)]))
  seqlevels(bamwhich, pruning.mode="coarse") <- seqlevels(BamFile(bamfile)) ## match seqlevels to bam file
  #return(bamwhich)
  #names(which) <- NULL
  #return(bamwhich)
  ## scanBam parameters
  bamparam.wins <-
    ScanBamParam(
      what = c("qwidth", "cigar"),
      ## qwidth is essential for FP size calculations
      flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE)
    )
  #return(bamparam.wins)
  ## determine if paired  
  is_paired <-
    suppressMessages(testPairedEndBam(BamFile(bamfile)))
  #return(is_paired)
  ## get GRanges
  if (is_paired) {
    gr <-
      GRanges(readGAlignmentPairs(
        file = bamfile,
        use.names = T,
        param = bamparam.wins
      ))
  } else{
    gr <-
      GRanges(readGAlignments(
        file = bamfile,
        use.names = T,
        param = bamparam.wins
      ))
  }
  ## update seqlevels
  seqlevels(gr, pruning.mode = "coarse") <-
    seqlevels(ex_by_tx)
  ### map read starts to transcripts
  tx <-
    suppressWarnings(mapToTranscripts(resize(gr, 1, "start"), ex_by_tx)) ## lift read starts to tx space
  
  ## update metadata
  mcols(tx) <-
    cbind(mcols(tx), mcols(gr[tx$xHits]))
  
  ## update strand
  strand(tx) <-
    "+"  ## set strand to + as mapToTranscripts retains genomic strand
  
  ## update seqlengths to get full coverage
  seqinfo(tx) <-
    Seqinfo(seqnames =  seqlevels(tx),
            seqlengths = txmodel$tx_len[match(seqlevels(tx), txmodel$tx_name)])
  
  ## test periodicity
  tx$frame <-
    get_frame(tx)
  per <-
    as.data.table(mcols(tx))[, .N, by = c("qwidth", "frame")]
  
  ## map to windows
  tx.win <-
    suppressWarnings(mapToTranscripts(tx, wins)) 
  mcols(tx.win) <-
    cbind(mcols(tx.win), mcols(tx[tx.win$xHits])[, -1][, -1]) ## transfer mcols but without xHits and transcriptHits
  
  # browser()
  
  ## determine position of maximum coverage in windows
  cov.dt <-
    cbind(data.table(pos = start(tx.win),
                     winclass = as.character(gsub("_.*$", "", seqnames(tx.win)))
    ),
    as.data.table(mcols(tx.win))
    )
  
  cov.dt[, pos := pos - (rwinsize +1) ] ## shift zero coordinate to center of wins
  
  ### calculate coverage
  cov.agg <-
    cov.dt[, .(cvg = .N), by = c("winclass", "pos", "qwidth")] ## calculate coverage
  
  ### determine offsets at start codons
  ribodat <-
    list(cov.agg = cov.agg,
         per = per)
  
  return(ribodat)
}