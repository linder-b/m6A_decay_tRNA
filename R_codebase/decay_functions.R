
#' A function to normalize time series expression data to t0.
#'
#' @param dt A data.table in long format with gene expression in gene_col and timepoints in time_col. Additional variables can be defined in id.vars. 
#' @param time_col Column name of time values. Should be numeric and the smallest value will be used as t0
#' @param expression_col Column name of gene expression values.
#' @param id.vars Gene expression is aggregated by these variables, so these should include genes, genotypes, treatments, etc.
#'   Note that if replicates exist and are not defined in id.vars, mean expression will be used for replicates, too.
#' @return A data table with the id.vars columns, the time_col and the expression_col (which is now normalized to t0).
norm_t0 <-
  function(dt = NULL, time_col = "time", expression_col = "expression", id.vars = c("gene", "treatment")) {
    if (!time_col %in% names(dt))
      stop (paste("Time column", time_col, "not found in table."))
    if (!expression_col %in% names(dt))
      stop (paste("Expression column", expression_col, "not found in table."))
    timepoints = unique(dt[[time_col]])
    tpt_char = as.character(timepoints)
    tpt_num = as.numeric(tpt_char)
    t0_num = tpt_num[which.min(tpt_num)]
    t0_char = tpt_char[which.min(tpt_num)]
    cat(paste("Normalizing to time", t0_char, "\n"))
    if (length(id.vars) > 0) {
      fml = paste(paste(id.vars, collapse = " + "), "~", time_col)
    } else{
      fml = paste(id.vars, "~", time_col)
    }
    
    dt <-
      dcast(dt, formula(fml), value.var = expression_col, fun.aggregate = mean)
    dt[, (tpt_char) := lapply(.SD, function(x) {
      x / get(t0_char)
    }), .SDcols = tpt_char]
    dt <-
      melt(
        dt,
        id.vars = id.vars,
        variable.factor = FALSE,
        variable.name = time_col,
        value.name = expression_col
      )
    if (is.numeric(timepoints))
      dt[, (time_col) := as.numeric(get(time_col))]
    return(dt[])
  }

#' A function to fit an exponential decay model.
#'
#' @param dt A data.table in long format with gene expression values and timepoints.
#' @param time_col Character string defining the column used as timepoints. Should be numeric and the smallest value will be used as t0
#' @param expression_col Character string defining the column with gene expression values.
#' @param method Either nls or lm.
#' @param verbose Whether errors in model fits are reported.
#' @return
#' A list containing:
#'   * the model fit
#'   * a table with the coefficients a and k and their standard errors
fit_decay <-
  function(dt = NULL,
           time_col = "time",
           expression_col = "expression",
           method = "nls",
           verbose = FALSE) {
    if (!time_col %in% names(dt))
      stop(paste("Time column", time_col, "not found in table!"))
    if (!expression_col %in% names(dt))
      stop(paste("Expression column", expression_col, "not found in table!"))
    if (!is.numeric(dt[[time_col]]))
      stop(paste(time_col, "must be numeric!"))
    if (!is.numeric(dt[[expression_col]]))
      stop(paste(expression_col, "must be numeric!"))
    
    if (method == "nls") {
      
      fml = paste(expression_col, "~ I(a * exp(-k *", time_col, "))")
      
      max_expr = max(dt[[expression_col]][is.finite(dt[[expression_col]]) &
                                            !is.na(dt[[expression_col]])])
      
      if (inherits(try(nls(
        formula =  fml,
        start = list(a = max_expr, k = 0),
        data = dt
      ),
      silent = !verbose)
      , "try-error")) {
        fit = NA
        a = as.numeric(NA)
        a_se = as.numeric(NA)
        k = as.numeric(NA)
        k_se = as.numeric(NA)
      } else {
        fit = nls(formula =  fml,
                  start = list(a = max_expr, k = 0),
                  data = dt)
        smr = summary(fit)
        a = smr$parameters["a", "Estimate"]
        a_se = smr$parameters["a", "Std. Error"]
        k = smr$parameters["k", "Estimate"]
        k_se = smr$parameters["k", "Std. Error"]
      }
    } else if(method == "lm") {
      
      #fml = paste(expression_col, "~ I(a -k *", time_col, "))") ## for use with nls
      fml = paste("log(", expression_col, ")",  "~", time_col)
      
      dt <-
        dt[!is.na(get(expression_col)) & is.finite(get(expression_col)), ]
      
      if (inherits(try(lm(
        formula =  fml,
        data = dt,
        na.action = na.exclude
      ),
      silent = !verbose)
      , "try-error")) {
        fit = NA
        a = as.numeric(NA)
        a_se = as.numeric(NA)
        k = as.numeric(NA)
        k_se = as.numeric(NA)
      } else {
        fit = lm(formula =  fml,
                 data = dt,
                 na.action = na.exclude)
        smr = summary(fit)
        a    = if (inherits(try(smr$coefficients["(Intercept)", "Estimate"], silent = !verbose), "try-error")){as.numeric(NA)}else{smr$coefficients["(Intercept)", "Estimate"]}  
        a_se = if (inherits(try(smr$coefficients["(Intercept)", "Std. Error"], silent = !verbose), "try-error")){as.numeric(NA)}else{smr$coefficients["(Intercept)", "Std. Error"]}
        k    = if (inherits(try(smr$coefficients[time_col, "Estimate"], silent = !verbose), "try-error")){as.numeric(NA)}else{smr$coefficients[time_col, "Estimate"]}
        k = -k ## we need to inverse the slope of lm
        k_se = if (inherits(try(smr$coefficients[time_col, "Std. Error"], silent = !verbose), "try-error")){as.numeric(NA)}else{smr$coefficients[time_col, "Std. Error"]} 
      }
    } else {
      stop("Undefined method!")
    }
    coef <-
      data.table(
        a = a,
        a_se = a_se,
        k = k,
        k_se = k_se)
    coef[, method := method]
    return(list(
      fit = fit,
      coef = coef))
  }

#' Predict expression decay from model coefficients.
#'
#' @param dt A data.table containing the coeffiecients, e.g. produced by fit_decay function. 
#' @param timepoints The timepoints to predict. 
#' @param time_col Timepoint column name.
#' @param expression_col Expression column name.
#' @return A data table with added columns for timepoints and predicted expression.
predict_decay <-
  function(dt = NULL, timepoints = seq(0, 6, 0.25), time_col = "time", expression_col = "expression") {
    idx = sort(rep(1:nrow(dt), length(timepoints))) ## expand table rows for all timepoints
    dt2 = dt[idx, ]
    dt2[, (time_col) := rep(timepoints, nrow(dt))]
    dt2[, (expression_col) := ifelse(method == "lm",
                                     exp(a - k * get(time_col)),
                                     ifelse(method == "nls", a * exp(-k * get(time_col)), NA))]
    return(dt2[])
  }


## eq. 13 from Jürges et al. 2018
decay_from_NTR <-
  function(label_time = NULL, NTR = NULL){
    dr =  -1*(1/label_time)*log(1-NTR)
    return(dr)
  }

t12_from_NTR <-
  function(label_time = NULL, NTR = NULL){
    dr =  -1*(1/label_time)*log(1-NTR)
    t12 = log(2)/dr
    return(t12)
  }

## check with data from fig. 5a (3h labelling of RNA with half-life of 2h results in NTR of 0.646):
t12_from_NTR(label_time = 3, NTR = 0.646)

## eq. 14 from Jürges et al. 2018
ntr_from_t12 <-
  function(label_time = NULL, t12 = NULL){
    decay_rate = log(2)/t12
    ntr = 1-(exp(1)^(-label_time*decay_rate))
    return(ntr)
  }

ntr_from_decay <-
  function(label_time = NULL, decay_rate = NULL){
    ntr = 1-(exp(1)^(-label_time*decay_rate))
    return(ntr)
  }

## common conversions
decay_from_t12 <-
  function(t12 = NULL){
    dr =  log(2)/t12
    return(dr)
  }

t12_from_decay <-
  function(dr = NULL){
    t12 =  log(2)/dr
    return(t12)
  }


## check with data from fig. 5a (3h labelling of RNA with half-life of 2h results in NTR of 0.646):
#ntr_from_t12(label_time = 3, t12 = 2)

#' A function to calculate decay fits.
#'
#' @param dt A data.table in long format with gene expression values. Timepoints, 
#' @param gene_col Character string defining the column used as gene identifiers
#' @param time_col Character string defining the column used as time points. Should be numeric and the smallest value will be used as t0
#' @return A list with two elements: The normalized expression relative to t0 and a data.table with the slopes (and the model fits)
#' @examples
get_decay_rates <-
  function(dt = NULL,
           gene_col = "gene",
           time_col = "timepoint",
           expression_col = "value",
           id.vars = NULL) {
    if (!is(dt, "data.table"))
      stop("data.table input required")
    missingCols = setdiff(c(gene_col, time_col, expression_col), names(dt))
    if (length(missingCols > 0))
      stop(paste(missingCols, "are required but not found in ", dt))
    ## copy to avoid back-propagation to dt outside of function
    dt <- copy(dt)
    if("timepointNames" %in% names(dt))stop('Column name "timepointNames" is reserved for internal use!')
    ## to use as column names in the wide table, we convert time to a character string (required for time with decimal points)
    dt[, timepointNames := paste0("time", as.character(get(time_col)))]
    timepoints <-
      gtools::mixedsort(unique(dt$timepointNames))
    #timepoints <- paste0("time", as.character(sort(unique(dt[[time_col]]))))
    
    
    if (length(id.vars) > 0) {
      fml = paste(paste(c(gene_col, id.vars), collapse = " + "), "~", "timepointNames")
    } else{
      fml = paste(gene_col, "~", "timepointNames")
    }
    # return(fml)
    
    ## check for duplicates (e.g. due to missing id.vars or duplicate gene identifiers)
    dups = dt[, .N, by = c(gene_col, time_col, expression_col, id.vars) ][ (N > 1), ]
    if(nrow(dups) > 0){warning("Duplicate entrys detected!"); return(dups)}
    
    ## normalize to first timepoint (note that replicate variance is lost when replicates are in id.vars)
    ### cast to wide for timepoints
    dt.norm <-
      dcast(data = dt,
            formula = as.formula(fml),
            value.var = expression_col)
    ### replace 0s in expression with the minimal non-zero (assuming this is detection limit) to enable log transformation for genes that drop below det. limit
    dt.norm[, timepoints[-1] := lapply(.SD, function(x){ifelse(x==0, min(x[x !=0 & !is.na(x)]), x)}), .SDcols = timepoints[-1]]
    ### normalize to first timepoint
    dt.norm[, timepoints[] := lapply(.SD, function(x){x / get(timepoints[1])}), .SDcols = timepoints[] ] ## square brackets are required for timepoints, as it is otherwise interpreted as column?!?
    ### melt to long
    dt.norm <-
      melt(dt.norm, id.vars = c(gene_col, id.vars), variable.factor = FALSE, variable.name = "timepointNames", value.name = "expression_rel_t0")
    ### make time numeric for model fit
    dt.norm[, time := as.numeric(gsub("^time", "", timepointNames))]
    
    ## get slopes
    slopes <-
      dt.norm[is.finite(log(expression_rel_t0)), {
        fit = lm(log(expression_rel_t0) - log(1) ~ 0 + time ) ## as we have normalized to t0, we set the intercept to 1 (cave: log) from: https://stackoverflow.com/a/7333292
        slope = coef(fit)[["time"]]
        list(slope = slope, model = list(fit) )
      },  by = c(gene_col, id.vars)]
    return(list(norm_t0 = dt.norm, slopes = slopes))
  }



## a function to correct for negative decay rates: we shift the decay rates so most of them are above 0;
## to take care of extreme outliers that would lead to unreasonable shifts, we use the 99% quantile
shift_ntile_zero <-
  function(x = NULL, ntile = 0.01){
    if(quantile(x, ntile, na.rm = TRUE) < 0 ){
      return( x + (-1 * quantile(x, ntile, na.rm = TRUE)) )
    } else {
      return(x)
    }
  }


#' A function to calculate the Codon Stabilization Coefficient (CSC).
#'
#' @param decayRates A data.table with decay rate data, 1st column is gene IDs, all others are considered decay rates (column names are used as dataset IDs).
#' @param codonUsage A data.table with codon usage data, 1st column is gene IDs, all others are considered codon usage (column names are used as codon IDs).
#' @param useCodonFrequency Whether codon counts or codon frequencies are used. Note that redundant codon labels will lead to double counting!
#' @param cscCorrelation Correlation method to use [spearman|pearson].
#' @param codonZeroFilter Whether non-existent codons in a gene should be removed before testing correlation.
#' @param ignoreNegativeDecay 
#' @param useDecayRates 
#' @param useRegression Whether to use the slope of a linear regression instead of the correlation coefficient. [TRUE|FALSE]
#' @return 
#' @examples
getCSC <-
  function(decayRates = NULL,
           codonUsage = NULL,
           useCodonFrequency = FALSE,
           cscCorrelation = "spearman",
           codonZeroFilter = TRUE,
           ignoreNegativeDecay = FALSE,
           useDecayRates = FALSE,
           useRegression = FALSE){
    
    idColName = names(decayRates)[1]
    
    if(idColName != names(codonUsage)[1] ){stop("Gene ID columns do not have the same name.")}
    commonIDs = intersect(decayRates[[1]], codonUsage[[1]])
    uniqueIDs = unique(c(decayRates[[1]], codonUsage[[1]]))
    
    if(length(commonIDs) < 1000 ){stop("Less than 1000 common gene IDs. Did you use the same IDs?")}
    
    ## cast to long
    gcu <-
      melt(codonUsage, id.vars = idColName, variable.name = "codonName", value.name = "codonOccurrence")
    ## apply filter
    if (codonZeroFilter == TRUE) {
      gcu = gcu[codonOccurrence > 0, ]
    } else if (codonZeroFilter == FALSE) {
      
    } else {
      stop("codonZeroFilter not defined (use [TRUE|FALSE]).")
    }
    
    ## calculate codon Frequency
    if (useCodonFrequency == TRUE) {
      ### get total number of codons in gene
      ### !!! beware of redundant codon labels !!!  (e.g. codons labelled with "all" and the DRACH classes)
      ### this will double-count codons
      ### set codonOccurrence to numeric, as it might be an integer and this throws a warning for each gene in the next step
      gcu[, codonOccurrence := as.numeric(codonOccurrence)]
      ### normalize codon occurrence by total number of codons in gene
      gcu[, codonOccurrence := codonOccurrence / sum(codonOccurrence),  by = idColName]
    } else if (useCodonFrequency == FALSE) {
      
    } else {
      stop("useCodonFrequency not defined (use [TRUE|FALSE]).")
    }
    
    
    gcu <-
      merge(gcu, decayRates, idColName)
    
    #return(gcu)
    
    gcu <-
      melt(gcu, id.vars = c(idColName, "codonName", "codonOccurrence"), variable.name = "decayDataName", value.name = "decayRate")
    
    if(ignoreNegativeDecay){gcu[decayRate <= 0, decayRate := NA]}
    
    if(useDecayRates == TRUE){
      ## use decay rate; to remain consistent with the half-life-based CSC where lower CSC is associated with less optimal / unstable, we use negative the decay rate
      gcu[, decayMetric := - decayRate ]
    } else if (useDecayRates == FALSE) {
      ## use half-life
      gcu[, decayMetric := log(2)/decayRate]
    } else {stop("Please define whether to use decay rate or half-life")}
    # return(gcu)
    
    ## calculate CSC
    
    if(useRegression == TRUE){
      csctab <-
        gcu[!is.na(decayMetric) & is.finite(decayMetric), .("csc" = coef(lm(decayMetric ~ codonOccurrence, .SD))[["codonOccurrence"]] ), by = c("codonName", "decayDataName") ]
    } else {
      csctab <-
        gcu[, .("csc" = suppressWarnings(cor.test(codonOccurrence, decayMetric, method = cscCorrelation)$estimate)), by = c("codonName", "decayDataName") ]
    }
    
    
    
    csctab <-
      dcast(csctab, codonName ~ decayDataName, value.var = "csc")
    
    return(csctab)
  }


read_grandslam_introns <-
  function(grandslam_tsv = NULL){
    dt <-
      fread(grandslam_tsv)
    setnames(dt, apply(str_split_fixed(basename(names(dt)), " ", 2), 1, function(x)paste(x, collapse= "##")))
    setnames(dt, gsub("##$", "", names(dt)) )
    dt[, "Symbol" := NULL]
    dt[, region := ifelse(grepl("_intronic", Gene), "intronic", "exonic")]
    dt[, gene_id := gsub("_intronic$", "", Gene)][, "Gene" := NULL]
    dt[, length := as.numeric(ifelse(Length == "null", NA, Length)) ][, "Length" := NULL]
    setcolorder(dt, c("gene_id", "length", "region"))
    dt <-
      melt(dt, id.vars = c("gene_id", "length", "region"), variable.factor = FALSE)
    dt[, run := gsub("##..*$", "", variable)]
    dt[, slam_var := gsub("^..*##", "", variable)]
    dt <-
      dcast(dt, gene_id + length + region + run ~ slam_var, value.var = "value")
    
    dt[, NTR := MAP]
    dt[, RPM := Readcount/sum(Readcount)*1e6, by = "run"]
    dt[region == "exonic", TPM  := countToTpm(counts = Readcount, effLen = length), by = "run" ]
    dt[, conv_rate := Conversions/Coverage ]
    
    return(dt)
  }
