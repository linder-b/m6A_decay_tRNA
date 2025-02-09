## extract bcr barcode info see: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
get_barcode_table <-
  function(bcr_full_barcode, include_full_barcode = TRUE) {
    require(data.table)
    require(tidyr)
    barcode_table <-
      as.data.table(do.call(rbind, strsplit(x = bcr_full_barcode, split = "-")))[, 1:7]
    setnames(
      barcode_table,
      c(
        "project",
        "tss",
        "participant",
        "sample_vial",
        "portion_analyte",
        "plate",
        "center"
      )
    )
    barcode_table[, sample := as.integer(gsub("[A-Z]$", "", sample_vial))]
    barcode_table[, bcr_sample_barcode := paste(project, tss, participant, sample, sep="-") ]
    barcode_table[, sample_type := ifelse(sample %in% 0:9,
                                          "Tumor",
                                          ifelse(
                                            sample %in% 10:19,
                                            "Normal",
                                            ifelse(sample %in% 20:29, "control", NA)
                                          ))]
    
    barcode_table[, submitter_id := paste(project, tss, participant, sep="-") ] ## submitter_id is used in TCGA clinical files
    if(include_full_barcode == TRUE){
      barcode_table$bcr_full_barcode <- bcr_full_barcode
      setcolorder(barcode_table, c("bcr_full_barcode", setdiff(names(barcode_table), "bcr_full_barcode")))
    }
    return(barcode_table)
  }

get_bcr_sample_type <-
  function(bcr_full_barcode, type = "Short.Letter.Code") {
    require(TCGAutils)
    data("sampleTypes")
    sample_key = gsub("[A-Z]$", "", do.call(rbind, strsplit(x = bcr_full_barcode, split = "-"))[, 4])
    sample_type = setkey(as.data.table(sampleTypes), Code)[sample_key, ][[type]]
    return(sample_type)
  }


get_bcr_patient_barcode <-
  function(bcr_full_barcode = NULL){
    #require(gsubfn)
    ## regex from:https://stackoverflow.com/a/28641543
    # bcr_patient_barcode = strapplyc(bcr_full_barcode, pattern = "^(?:[^-]*-){2}(?:[^-]*){1}", simplify = c)
    bcr_patient_barcode = apply(matrix(do.call(rbind, strsplit(bcr_full_barcode, split = "-"))[, 1:3], ncol = 3), 1, function(x){paste(x, collapse = "-")} )
    return(bcr_patient_barcode)
  }

get_bcr_sample_barcode <-
  function(bcr_full_barcode = NULL){
    # require(gsubfn)
    ## regex from:https://stackoverflow.com/a/28641543
    # bcr_sample_barcode = gsub("[A-Z]$", "", strapplyc(bcr_full_barcode, pattern = "^(?:[^-]*-){3}(?:[^-]*){1}", simplify = c))
    bcr_sample_barcode = apply(do.call(rbind, strsplit(bcr_full_barcode, split = "-"))[, 1:4], 1, function(x){gsub("[A-Z]$", "", paste(x, collapse = "-"))} )
    return(bcr_sample_barcode)
  }



get_tcga_clinical_table <-
  function(project = NULL, tmpDir = NULL) {
    tcgaTableFile = file.path(tmpDir, paste0("TCGAquery.", project, ".clinical.RData"))
    if (!file.exists(tcgaTableFile)) {
      data <-
        as.data.table(GDCquery_clinic(project = project, type = "clinical"))[, project_id := project]
      saveRDS(data, tcgaTableFile)
    } else{
      data = readRDS(tcgaTableFile)
    }
    return(data)
  }


## a wrapper to download GDC query data
get_gdc_query <-
  function(query_params = NULL,
           destdir = NULL) {
    tableFile =  file.path(tmpDir, paste0(digest::digest(query_params), ".table.RData"))
    
    if (!file.exists(tableFile)) {
      ## generate query
      queryFile = file.path(tmpDir, paste0(digest::digest(query_params), ".GDCquery.RData"))
      if (!file.exists(queryFile)) {
        query <- do.call(GDCquery, query_params)
        saveRDS(object = query, file = queryFile)
      } else{
        query <- readRDS(queryFile)
      }
      
      ## download data
      files <-
        with(query$results[[1]],
             gsub(
               " ",
               "_",
               paste(
                 tmpDir,
                 project,
                 ifelse(query$legacy, "legacy", "harmonized"),
                 data_category,
                 data_type,
                 file_id,
                 file_name,
                 sep = "/"
               )
             ))
      if (!all(file.exists(files))) {
        GDCdownload(query = query,
                    directory = tmpDir)
      }
      
      ## generate table
      data <-
        foreach(f = files, .combine = cbind) %do% {
          fread(f, select = 2)
        }
      data <-
        cbind(fread(files[1], select = 1), data)
      setnames(data, c("ensembl_gene_id", query$results[[1]]$cases))
      saveRDS(data, tableFile)
    } else{
      data <- readRDS(tableFile)
    }
    return(data)
  }


get_recount <-
  function(study_id = NULL, temp_dir = NULL){
    require(recount)
    require(SummarizedExperiment)
    if(!dir.exists(temp_dir)){dir.create(tmp_dir, recursive = TRUE)}
    study_file = file.path(temp_dir, file.path(study_id, "rse_gene.Rdata"))
    if (!file.exists(study_file)) {
      download_study(project = study_id,
                     type = "rse-gene",
                     outdir = dirname(study_file))
    }
    load(study_file) ## loads rse_gene object
    
    rse_gene <-
      scale_counts(
        rse = rse_gene,
        by = "mapped_reads",
        targetSize = 4e+07,
        L = 100,
        factor_only = FALSE,
        round = TRUE
      )
    
    return(rse_gene)  
    
  }