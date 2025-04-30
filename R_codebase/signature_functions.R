require(survival)



get_best_correlated <-
  function(gene, data) {
    corrs = cor(data)[, gene]
    corrs = corrs[!grepl(gene, names(corrs))]
    return(corrs[which.max(abs(corrs))])
  }


get_signature_model <-
  function(genes = NULL,
           data = NULL,
           fml.lhs = NULL,
           ncores = 1) {
    var_lhs = unlist(strsplit(gsub("Surv\\(|\\)", "", fml.lhs), split = ", "))
    data = data[, c(var_lhs, genes), with = FALSE]
    registerDoParallel(cores = ncores)
    coxUniList <-
      foreach(gn = genes) %dopar% {
        coxph(formula(paste(fml.lhs, "~", gn)),
              data = data,
              x = TRUE)
      }
    registerDoSEQ()
    coxUniTab <-
      data.table(
        gene = genes,
        coef = unlist(lapply(coxUniList, function(x)
          coef(x)[[1]])),
        waldtest = unlist(lapply(coxUniList, function(x)
          summary(x)[["waldtest"]][["pvalue"]])),
        logtest  = unlist(lapply(coxUniList, function(x)
          summary(x)[["logtest"]][["pvalue"]]))
      )
    setkeyv(coxUniTab, "gene")
    topGenes = coxUniTab[waldtest < 0.05, gene]
    #topGenes = coxUniTab[["gene"]][order(coxUniTab[["waldtest"]])]
    coxMulti <-
      coxph(formula = formula(paste(
        fml.lhs, "~", paste(topGenes, collapse = " + ")
      )),
      data = data,
      x = TRUE)
    results = setnames(as.data.table(coef(summary(coxMulti)), keep.rownames = TRUE), "rn", "feature")[0, ]
    coxIterationList = list(starting = coxMulti)
    for(gn in topGenes){
      cat(paste0("processing ", gn, ": "))
      if(gn %in% names(coef(coxMulti))){
        if(coef(summary(coxMulti))[gn, 'Pr(>|z|)'] < 0.05){
          cat(paste("is significant!", "\n"))
          results = rbind(results,
                          setnames(as.data.table(coef(summary(coxMulti)), keep.rownames = TRUE), "rn", "feature")[ feature == gn, ])
        } else {
          best = get_best_correlated(gene = gn, data = data[, names(coef(coxMulti)), with = FALSE])
          new_feature = paste0(gn, "_", names(best))
          cat(paste("not significant; merging", gn, "and", names(best), "(cor=", round(best, 3), ")","\n"))
          data[, eval(new_feature) := rowMeans(as.matrix(.SD)), .SDcols = c(gn, names(best))]
          #data[[new_feature]] = ifelse(best >0, data[[gn]] + data[[names(best)]], data[[gn]] - data[[names(best)]] )
          coxMulti <-
            coxph(formula = formula(paste(
              fml.lhs, "~", paste( c(new_feature, setdiff(names(coef(coxMulti)), c(gn, names(best) )) ) , collapse = " + ")
            )),
            data = data,
            x = TRUE)
          coxIterationList[[gn]] = coxMulti
          results = rbind(results,
                          setnames(as.data.table(coef(summary(coxMulti)), keep.rownames = TRUE), "rn", "feature")[ feature == new_feature, ])
        }
      }else{cat(paste("has been removed in prior iteration", "\n"))}
    }
    
    return(list(coxUniTab  = coxUniTab,
                coxIter   = coxMulti,
                coxData    = data[, setdiff(names(coef(coxMulti)), genes), with = FALSE],
                coxResults = results,
                coxIterationList = coxIterationList,
                coxUniList = coxUniList))
  }


## calculate signature from existing iterative survival model and new expression table dt (genes in rows, samples in columns) using the coefficients of the iterative model
getIterSignature <-
  function(signature_model = NULL,
           dt = dt) {
    ## extract list of the new features (of combined genes) used in the iterative model
    newFeatureList <-
      strsplit(names(coef(signature_model$coxIter)), split = "_")
    names(newFeatureList) <-
      names(coef(signature_model$coxIter))
    ## adapt names to current data (mouse)
    newFeatureList <-
      get_geneset(gsList = newFeatureList, dict = unique(dt$gene))
    ## feature expression table (median expression as used in model)
    featureExprData <-
      cbind(sample = names(dt[,-1]),
            foreach(
              nfn = names(newFeatureList),
              nf = newFeatureList,
              .combine = cbind
            ) %do% {
              setnames(data.table(colMedians(as.matrix(dt[gene %in% nf,-1]), na.rm = TRUE)), "V1", nfn)
            })
    signature <-
      Reduce('+', lapply(names(newFeatureList), function(x) {
        featureExprData[[x]] *  coef(signature_model$coxIter)[x]
      }))
    
    return(data.table(sample = featureExprData$sample, signatureIter = signature ))
  }

getMultiSignature <-
  function(signature_model = NULL,
           dt = dt) {
    ## extract list of the new features (of combined genes) used in the iterative model
    newFeatureList <-
      strsplit(names(coef(signature_model$coxIterationList$starting)), split = "_")
    names(newFeatureList) <-
      names(coef(signature_model$coxIterationList$starting))
    ## adapt names to current data (mouse)
    newFeatureList <-
      get_geneset(gsList = newFeatureList, dict = unique(dt$gene))
    ## feature expression table (median expression as used in model)
    featureExprData <-
      cbind(sample = names(dt[,-1]),
            foreach(
              nfn = names(newFeatureList),
              nf = newFeatureList,
              .combine = cbind
            ) %do% {
              setnames(data.table(colMedians(as.matrix(dt[gene %in% nf,-1]), na.rm = TRUE)), "V1", nfn)
            })
    signature <-
      Reduce('+', lapply(names(newFeatureList), function(x) {
        featureExprData[[x]] *  coef(signature_model$coxIterationList$starting)[x]
      }))
    
    return(data.table(sample = featureExprData$sample, signatureMulti = signature ))
  }


getUniSignature <-
  function(signature_model = NULL,
           dt = dt) {
    ## adapt gene names to current data (mouse), then unlist (as we are directly using gene symbols)
    coxUniGenes <-
      # translate_symbol(signature_model$coxUniTab$gene, dict = unique(dt$gene), fill = TRUE)
    #get_geneset(list("coxUniGenes" = signature_model$coxUniTab$gene))
    unlist(get_geneset(setNames(as.list(signature_model$coxUniTab$gene), signature_model$coxUniTab$gene), dict = unique(dt$gene) ))
    ## feature expression table (just transpose dt as we are directly using gene symbols)
    featureExprData <-
      setnames(tdt(dt), 1, "sample")
    signature <-
      Reduce('+', lapply(names(coxUniGenes), function(x) {
        featureExprData[[coxUniGenes[[x]] ]] *  setkeyv(signature_model$coxUniTab, "gene")[x, coef]
      }))
    return(data.table(sample = featureExprData$sample, signatureUni = signature ))
  }

#' Calculate gene signature.
#' 
#' @param dt Expression table. A data.table, first column is sample IDs, other columns are gene expression values (other column names are gene IDs).
#' @param sig_coef Coefficient table. A data.table, first column is gene IDs, second column is model coefficients.
#' @returns A data.table, first column is sample IDs, second column is gene signature.
#' @examples
#' get_signature(dt, sig_coef)
get_signature <-
  function(dt = NULL, sig_coef = NULL){
    if(length(names(dt)) != length(unique(names(dt))) ){
      stop("Duplicate gene names in sig_coef")
    }
    
    model_genes = intersect(sig_coef[[1]], names(dt))
    
    missing_genes = setdiff(sig_coef[[1]], names(dt))
    
    setkeyv(sig_coef, names(sig_coef[, 1]))
    
    if(length(missing_genes > 0)){
      warning(paste("The following", length(missing_genes), "genes are not found in the expression table and will be ignored:", paste(missing_genes, collapse = ", ") ))
    }
    dt = dt[, c(1, which(names(dt) %in% model_genes)), with = FALSE]
    
    dtm = melt(dt, id.vars = names(dt[, 1]), variable.factor = FALSE, variable.name= names(sig_coef[, 1]))
    
    dtm = merge(dtm, sig_coef)
    
    sig = dtm[, .(signature = sum(value*coef)), by = names(dt[, 1])]
    
    sig = setkeyv(sig, names(dt[, 1]))[dt[[1]], ]
    
    return(sig)
  }

cut_signature <-
  function(signature = NULL, g = 4, cuts = NULL){
    f <-
      cut2(x = signature, g = g)
    f <-
      ifelse(test = f == levels(f)[1],
             yes  = "low",
             no   = ifelse(test = f == levels(f)[length(levels(f))],
                           yes  = "high",
                           no   = "intermediate")
      )
    f <-
      factor(f, levels = c("low", "intermediate", "high"))
    f<-
      droplevels(factor(f))
    return(f)
  }

### a generalized function to calculate expression signatures using model coefficients and id conversions from tables with genes in columns
# getSignature <-
#   function(modelGenes = NULL,
#            modelCoefficients = NULL,
#            modelName = "signature",
#            convertIdFunction = NULL,
#            exprData = NULL,
#            return_table = FALSE) {
#     tbl = data.table(gn = modelGenes,
#                      cf = modelCoefficients,
#                      id = if(!is.null(convertIdFunction)){convertIdFunction(modelGenes)}else{modelGenes})
#     if(return_table){return(tbl)}
#     if(any(is.na(tbl$id))){warning(paste("No IDs found for", paste(tbl[is.na(id), gn], collapse = ", ") ))}
#     if( ! length(modelGenes) == nrow(tbl)){stop("Genes not represented by unique IDs")}
#     nomatch <-
#       setdiff(tbl$id, names(exprData))
#     if(length(nomatch) > 0){
#       tbl <-
#         tbl[! gn %in% nomatch, ]
#       warning(paste("The IDs", paste(nomatch, collapse = ", "), "were not found in the expression table."))
#     }
#     #return(tbl)
#     signature <-
#       Reduce('+', lapply(tbl[!is.na(id), id], function(x) {
#         exprData[[x]] *  setkeyv(tbl, "id")[x, cf]
#       }))
#     sigTab <-
#       data.table(sampleID = exprData[[1]], signature = signature)
#     setnames(sigTab, "sampleID", names(exprData[, 1]))
#     setnames(sigTab, "signature", modelName)
#     setkeyv(sigTab, names(exprData[, 1]))
#     return(sigTab$signature)
#   }

### a generalized function to calculate expression signatures using model coefficients and id conversions from tables with genes in columns
# getSignature <-
#   function(modelGenes = NULL,
#            modelGeneCoefficients = NULL,
#            modelName = NULL,
#            exprData = NULL) {
#     tbl = data.table(gn = modelGenes,
#                      cf = modelGeneCoefficients)
#     #return(tbl)
#     genes_not_found = setdiff(tbl$gn, names(exprData))
#     if(length(genes_not_found) > 0){
#       warning(paste("The following genes are not in the expression table:", paste(genes_not_found, collapse = ", "), "\n"))
#       tbl <-
#         tbl[! gn %in% genes_not_found, ]
#     }
#     signature <-
#       Reduce('+', lapply(tbl$gn, function(x) {
#         exprData[[x]] *  setkeyv(tbl, "gn")[x, cf]
#       }))
#     
#     sigTab <-
#       data.table(sampleID = exprData[[1]], signature = signature)
#     setnames(sigTab, "sampleID", names(exprData[, 1]))
#     setnames(sigTab, "signature", modelName)
#     setkeyv(sigTab, names(exprData[, 1]))
#     return(sigTab[, modelName, with = FALSE])
#   }


getSignature <-
  function(exprMat = NULL, signatureCoefficients = NULL) {
    signatureCoefficients <-
      as.data.table(signatureCoefficients)
    if(!"gene" %in% names(signatureCoefficients)){stop('"Requried column "gene" not found in signatureCoefficients table.')}
    if(!"coef" %in% names(signatureCoefficients)){stop('"Requried column "coef" not found in signatureCoefficients table.')}
    
    missingGenes <-
      setdiff(signatureCoefficients$gene, names(exprData))
    
    if(length(missingGenes) > 0){warning(paste("Genes", paste(missingGenes, collapse = ", "),"are not in the expression matrix!"))}
    
    sig <-
      apply(sweep(exprMat, 2, setkeyv(signatureCoefficients, "gene")[colnames(exprMat), coef], "*"), 1, sum )
    
    return(sig)
  }

calculate_mRNAsi_signature_coefficients <-
  function(tmpDir = getwd(), synapser_email = NULL, synapser_pwd = NULL, idtab = NULL ){
    
    
    ## from: http://tcgabiolinks.fmrp.usp.br/PanCanStem/mRNAsi.html
    
    # Malta et al. 2018. Machine Learning Identifies Stemness Features Associated with Oncogenic Dedifferentiation.
    # Cell. 173:338-354.e15. doi:10.1016/j.cell.2018.03.034.
    
    
    ## dependencies
    # deps <- c("gelnet","dplyr","gdata","DT")
    # for(pkg in deps)  if (!pkg %in% installed.packages()) install.packages(pkg, dependencies = TRUE)
    # library(devtools)
    # if (!"synapseClient" %in% installed.packages())
    #   install_github('Sage-Bionetworks/rSynapseClient', ref = 'develop')
    
    library(gelnet)
    library(dplyr)
    library(biomaRt)
    library(synapser)
    
    # Maps ENSEMBL IDs to HUGO
    # Use srcType = "ensembl_gene_id" for Ensembl IDs
    # Use srcType = "entrezgene" for Entrez IDs
    genes2hugo <- function( v, srcType = "ensembl_gene_id" )
    {
      ## Retrieve the EMSEMBL -> HUGO mapping
      ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl" )
      ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
      
      ## Make sure there was at least one mapping
      if( nrow(ID) < 1 ) top( "No IDs mapped successfully" )
      
      ## Drop empty duds
      j <- which( ID[,2] == "" )
      if( length(j) > 0 ) ID <- ID[-j,]
      stopifnot( all( ID[,1] %in% v ) )
      
      ID
    }
    
    # Load RNAseq data
    synapser::synLogin(email=synapser_email, password=synapser_pwd)
    synRNA <- synapser::synGet( entity = "syn2701943", downloadLocation = file.path(tmpDir))
    
    
    X <- read.delim( synRNA$path ) %>%
      tibble::column_to_rownames( "tracking_id" ) %>% as.matrix
    X[1:3,1:3]
    
    # Retrieve metadata
    synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
    
    Y <-
      fread(synMeta$filepath)
    Y[, UID := gsub("-|\\+", ".", UID)]
    
    Y[grepl("014BEB" , UID), ]
    Y[grepl("039ECTO", UID), ]
    
    # Retrieve the labels from the metadata
    #y <- Y[colnames(X),]
    y <- setkeyv(Y, "UID")[colnames(X), Diffname_short]
    
    names(y) <- colnames(X)
    # Fix the missing labels by hand
    y["SC11.014BEB.133.5.6.11"] <- "EB"
    y["SC12.039ECTO.420.436.92.16"] <- "ECTO"
    
    ## Drop the splice form ID from the gene names
    v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
    
    rownames(X) <- v
    
    # Map Ensembl IDs to HUGO
    #V <- genes2hugo( rownames(X) )
    V <-
      data.frame(setnames(setkeyv(unique(idtab[, c("ensembl_gene_id", "external_gene_name")]), "ensembl_gene_id")[rownames(X), ], "external_gene_name", "hgnc_symbol"))
    
    V <-
      V[complete.cases(V), ]
    
    X <- X[V[,1],]
    rownames(X) <- V[,2]
    X[1:3,1:3]
    
    fnGenes = NULL
    
    if(!is.null(fnGenes)){
      vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
      VE <- genes2hugo( vGenes, "entrezgene" )
      X <- X[intersect( rownames(X), VE[,2] ),]
    }
    
    
    
    m <- apply( X, 1, mean )
    X <- X - m
    X[1:3,1:3]
    
    j <- which( y == "SC" )
    X.tr <- X[,j]
    X.tr[1:3,1:3]
    
    X.bk <- X[,-j]
    X.bk[1:3,1:3]
    
    require(gelnet)
    
    system.time(
      mm <- gelnet( t(X.tr), NULL, 0, 1 )
    )
    
    return(mm$w)
    
  }

## because the paper supplement has only TCGA samples, we calculate mRNAsi for all samples
get_mRNAsi <-
  function(dt = NULL, fnSig = NULL) {
    ## dt is a data.table where rows are genes and columns are samples, except for first column, which has the gene IDs
    ## the model we load is trained in malta_2018_pancanstem.Rmd
    if(!file.exists(fnSig))stop("Signature file not found.")
    ## read in signature
    w <-
      read.delim(fnSig, header = FALSE, row.names = 1) %>% as.matrix() %>% drop()
    #w[1:10]
    
    ## matrix with genes in rows and samples in cols
    # X <-
    #   t(as.matrix(exprData[, -1], rownames = exprData[[1]]))
    X <-
      as.matrix(dt[, -1], rownames = dt[[1]])
    #X[1:3,1:3]
    
    ## if matrix has additional genes, remove them
    X <-
      X[(rownames(X) %in% names(w)), ]
    
    ## adapt signature to available genes
    stopifnot(all(rownames(X) %in% names(w)))
    w <- w[rownames(X)]
    #w[1:5]
    
    ## calc signature
    s <-
      apply(X, 2, function(z) {
        cor(z, w, method = "sp", use = "complete.obs")
      })
    #s[1:5]
    
    ## Scale the scores to be between 0 and 1
    s <- s - min(s)
    s <- s / max(s)
    #s[1:5]
    #return(s)
    return(setnames(as.data.table(s, keep.rownames = TRUE), c("sample", "mRNAsi")))
  }


get_proliferativeIndex <-
  function(dt = NULL) {
    ## dt is a data.table where rows are genes and columns are samples, except for first column, which has the gene IDs
    require(ProliferativeIndex)
    df = as.data.frame(as.matrix(dt[,-1], rownames = dt[[1]]))
    #ids = dt[[1]]
    #return(df)
    pin = calculatePI(readDataForPI(vstData = df, modelIDs = rownames(df)))
    setnames(as.data.table(pin, keep.rownames = TRUE), c("sample", "PI"))[]
  }

get_exbinFit <-
  function(dt = NULL, exbins = NULL){
    ## dt is a data.table where rows are genes and columns are samples, except for first column, which has the gene IDs
    ## exbins is a data table with gene IDs in first column and a column named "exbin" that has CDS exon length bins (e.g. exlenTabGene)
    if(! any(names(exbins) == names(dt[, 1])) ){stop(paste0('Column "', names(dt[, 1]), '" not found in exbins table.'))}
    if(! any(names(exbins) == "exbin") ){stop(paste("Column", "exbin", "not found in exbins table."))}
    merge(x = melt(data = dt,
                   id.vars = c(names(dt[, 1]) ),
                   variable.name = "sample",
                   variable.factor = FALSE, 
                   value.name = "expression"
    ),
    y = exbins[, c(names(dt[, 1]), "exbin"), with = FALSE],
    by = names(dt[, 1]))[, {
      
      fit = lm(expression ~ as.integer(exbin))
      coeff = summary(fit)$coefficients
      data.table(
        intercept = coeff["(Intercept)", "Estimate"],
        intercept_stderr = coeff["(Intercept)", "Std. Error"],
        intercept_pval = coeff["(Intercept)", "Pr(>|t|)"],
        slope = coeff["as.integer(exbin)", "Estimate"],
        slope_stderr = coeff["as.integer(exbin)", "Std. Error"],
        slope_pval = coeff["as.integer(exbin)", "Pr(>|t|)"]
      )
    }, by = "sample"]
  }

#### ----- plot functions -----

palval_signature =  c(low = "#FFDE0D", intermediate = "grey", high = "#6E94CD")

## colors as in Malta et al. Cell 2018
palval_pam50 = c("Normal tissue" = "grey90",
                 "LumA" = "#0039FF",
                 "LumB" = "#00CFFF",
                 "Normal-like" = "cyan",
                 "Her2" = "#FFBAC4",
                 "Basal" = "#FF0000", 
                 "claudin-low" = "orange")

## order from Ceccarelli et al. Cell 2016 Table 2
palval_gbmSubtype = c(
  "Codel" = "#FFFA00",
  "G-CIMP-low" = "#A2E1EC",
  "G-CIMP-high" = "#0039FF",
  "Classic-like" = "#FFAE00",
  "Mesenchymal-like" = "#C93FF8",
  "LGm6-GBM" = "#FF0000",
  "PA-like" = "#007600"
)

palval_gbmGrade = c(
  "G2" = "#00FD00",
  "G3" = "#FF0000",
  "G4" = "#0039FF"
)


# palval_sets = c(m6A     = "#FFDE0D",
#                 m6A_neg = "grey80",
#                 mcm5    = "grey60",
#                 s2U     = "#6E94CD",
#                 proliferation = "#F8766D")

signaturePlot <-
  list(
    geom_bar(stat="identity", col = NA, width = 1 ),
    scale_x_continuous(expand = c(0,0)),
    scale_fill_manual(values = palval_signature),
    #geom_hline(yintercept = 0),
    labs(x="Signature rank", y = "Signature", fill = "Bin"),
    theme_pubr(),
    theme(legend.pos = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()),
    NULL
  )

stripPlot <-
  list(
    geom_tile(),
    scale_x_continuous(expand = c(0,0)),
    scale_y_discrete(expand = c(0,0)),
    labs(x="Signature rank"),
    theme_pubr() +
      theme(legend.pos = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()),
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 12)),
    NULL
  )

## a function for combined plotting of rank and strips
# plot_rankstrip <-
#   function(survData = NULL, sigBinsN = 4) {
#     p.rank <-
#       ggplot(survData, aes(x= rank(signature, ties = "first"), y = signature, fill = cut_signature(signature = signature, g = sigBinsN) )) +
#       signaturePlot
#     p.pam50 <-
#       ggplot(survData, aes(x = rank(signature), y = "PAM50", fill = pam50)) +
#       stripPlot +
#       scale_fill_manual(values = palval_pam50)
#     # p.er <-
#     #   ggplot(survData, aes(x = rank(signature), y = "ER status", fill = er)) +
#     #     stripPlot +
#     #     scale_fill_manual(values = palval_ER)
#     p.all <-
#       ggarrange(p.rank +
#                   theme(legend.pos = "none", axis.title.x = element_blank()),
#                 p.pam50 +
#                   theme(legend.pos = "none", axis.title.x = element_blank()),
#                 # p.er +
#                 #   theme(legend.pos = "none", axis.title.x = element_blank()),
#                 align = "v", ncol = 1, heights = c(4, 1, 1)
#       )
#     return(p.all)
#   }

plot_rankstrip <-
  function(survData = NULL, signatureVar = "signature", signatureBinN = 4, groupVar = NULL, groupPalette = NULL) {
    tbl <-
      data.table(sig = survData[[signatureVar]])
    tbl[, rnk := rank(sig, ties = "first")]
    tbl[, bin := cut_signature(signature = sig, g = signatureBinN)]
    tbl <-
      cbind(tbl, survData[, groupVar, with = FALSE])
    #return(tbl)
    
    p.rank <-
      ggplot(tbl, aes(x= rnk, y = sig, fill = bin )) +
      signaturePlot +
      theme(legend.pos = "none", axis.title.x = element_blank())
    
    p.strip.list <-
      foreach(v = groupVar) %do% {
        ggplot(tbl, aes_string(x = "rnk", y = 1, fill = v), ) +
          stripPlot +
          theme(legend.pos = "none", axis.title.x = element_blank()) +
          #scale_y_discrete(expand = c(0,0), breaks = 1, labels = v, drop = FALSE) +
          if(!is.null(groupPalette)){scale_fill_manual(values = groupPalette)} else {NULL}
      }
      
    p.list <-
      c(list(p.rank),
        p.strip.list)
    
    p.all <-
      ggarrange(plotlist = p.list,
                align = "v", ncol = 1, heights = c(4, rep(1, length(p.strip.list)))
      )
    return(p.all)
  }

plot_rank <-
  function(survData = NULL, signatureVar = "signature", signatureBinN = 4, omit_legend = FALSE) {
    tbl <-
      data.table(sig = survData[[signatureVar]])
    tbl[, rnk := rank(sig, ties = "first")]
    tbl[, bin := cut_signature(signature = sig, g = signatureBinN)]
    #return(tbl)
    
    p <-
      ggplot(tbl, aes(x= rnk, y = sig, fill = bin )) +
      signaturePlot +
      theme(axis.title.x = element_blank())
    if(omit_legend == TRUE){
      p <- p + theme(legend.pos = "none")
    }
    return(p)
  }

plot_strip <-
  function(survData = NULL, signatureVar = "signature", signatureBinN = 4, groupVar = NULL, groupPalette = NULL) {
    tbl <-
      data.table(sig = survData[[signatureVar]])
    tbl[, rnk := rank(sig, ties = "first")]
    tbl[, bin := cut_signature(signature = sig, g = signatureBinN)]
    tbl <-
      cbind(tbl, survData[, groupVar, with = FALSE])
    #return(tbl)
    
    p <-
        ggplot(tbl, aes_string(x = "rnk", y = groupVar, fill = groupVar), ) +
          stripPlot +
          theme(legend.pos = "none", axis.title.x = element_blank()) +
          #scale_y_discrete(expand = c(0,0), breaks = 1, labels = groupVar, drop = FALSE) +
          if(!is.null(groupPalette)){scale_fill_manual(values = groupPalette)} else {NULL}
      
    return(p)
  }


## a function for KM plot
plot_km <-
  function(survData = NULL, sigBinsN = 4) {
    ggsurvplot(
      survfit(formula = formula(paste(fml.lhs, "~", "cut_signature(signature, g = sigBinsN)")),
              data = survData),
      data = survData,
      palette = as.character(palval_signature[levels(droplevels(cut_signature(survData$signature)))]),
      pval = TRUE,
      legend = "none",
      #xscale = "d_y"
    )$plot
  }  


pam50Violin <-
  function(ref.group = ".all."){
    list(
      geom_violin(fill = NA),
      geom_boxplot(width = 0.2, outlier.shape = NA),
      geom_quasirandom(alpha = 0.2),
      stat_summary(fun.y = mean, geom = "hline", aes(x = 1, yintercept = ..y..), linetype = 2, col = "blue"),
      stat_summary(fun.y = mean, geom = "point", col ="black", size = 2),
      stat_compare_means(ref.group = ref.group, label ="p.signif"),
      scale_color_manual(values = palval_pam50),
      #labs(title = "Signature", x = "", y = "Signature"),
      theme_pubr(legend = "none"),
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()),
      NULL
    )
  }
    

summarizeGenesetExpr <-
  ## takes an expression table (samples in rows, genes in cols, sampleIDs in first col) and a list of geneset vectors to calculate means or medians
  function(exprData = NULL,
           geneset = NULL,
           comb_method = "mean") {
    exprSetList <-
      foreach(gs = geneset, gsn = names(geneset)) %do% {
        colSelect = intersect(gs, names(exprData))
        setnames(data.table("V1" = exprData[, if (comb_method == "median") {
          rowMedians(as.matrix(.SD), na.rm = TRUE)
        } else if (comb_method == "mean") {
          rowMeans(.SD, na.rm = TRUE)
        } else if (comb_method == "sum") {
          rowSums(.SD, na.rm = TRUE)
        } else{
          stop("comb_method needs to be one on [mean|median|sum]")
        }, .SDcols = colSelect]), "V1", gsn)
      }
    
    exprSetData <-
      do.call(cbind, c(exprData[, 1], exprSetList))
    return(exprSetData)
  }

