# install ImmuneSpaceR and ImmuneSignatures from GitHub
# devtools::install_github("rglab/ImmuneSignatures", build_vignettes = TRUE)
# devtools::install_github("rglab/ImmuneSpaceR")

library(ImmuneSignatures) # to load data and helper functions
library(Biobase) # biocLite
library(preprocessCore) # biocLite
library(stringr)


# function prepare ExpressionSet object from ImmuneSpace data
prepEset <- function(em, adjHai, d0 = TRUE){
  
  # get geneSyms now because prep work requires them to be absent
  new_fData <- data.frame(gene_symbol = em$geneSymbol, stringsAsFactors = F)
  rownames(new_fData) <- rownames(em)
  
  # rm GeneSyms and unwanted timepoints
  em <- if( d0 ){ em[ , grep("d0", colnames(em)) ] }else{ em[, !(colnames(em) == "geneSymbol")] }
  
  # determine which em subs have hai data and subset
  haiSubs <- sapply(colnames(em), FUN = function(name){
    sub <- strsplit(name, split = "_", fixed = T)[[1]][1]
    return( sub %in% adjHai$subject )
  })
  
  em <- em[, haiSubs ]
  
  # if all timepoints, create new hai that has one row per timepoint
  if( !d0 ){
    adjHai <- rbindlist(lapply(colnames(em), FUN = function(name){
      tmp <- strsplit(name, split = "_", fixed = T)[[1]]
      row <- adjHai[ adjHai$subject == tmp[1], ]
      day <- gsub("neg", "-", tmp[2])
      row$timepoint.day <- gsub("d", "", day)
      row$subject <- name
      return(row)
    }))
  }
  
  # Get Subjects from em and ensure hai / em contain same subs in same order
  emSubs <- if( d0 ){
    str_match(colnames(em), "(((SUB)|(sub))_?[0-9]+)(\\.[0-9])?_(.*)")[, 2L]
  } else {
    colnames(em)
  }
  
  adjHai <- adjHai[ adjHai$subject %in% emSubs, ]
  adjHai <- adjHai[ order(match(adjHai$subject, emSubs)), ]
  rownames(adjHai) <- colnames(em)
  
  # Make eset
  eset <- ExpressionSet(assayData = as.matrix(em),
                        featureData = AnnotatedDataFrame(new_fData))
  pData(eset) <- droplevels(adjHai)
  
  return(eset)
}


studies <- c("SDY212","SDY404","SDY400")
eset_list <- list()
# sdy80all <- "" # need eset with all time points for some figures
dir.create("generated_data/HIPC")

for(sdy in studies){
  labkey.url.path <- paste0("/Studies/", sdy)
  con <- CreateConnection()
  mats <- con$data_cache$GE_matrices$name
  
  # --------GET GENE EXPRESSION DATA------------
  
  # "ImmSig" annotation is particular to time the assays were read and developed
  # from the original manuscript's files, which were in turn often generated during
  # the sequencing process. Using the default / BioConductor Anno from ImmuneSpace
  # causes slightly different results and therefore is not used.
  
  # OutputType must be "raw" in order to pull non-normalized data for Studies
  # 212, 63, 404, and 67.  SDY80 and SDY400 did not have non-normalized available
  # from GEO with sampleNames that could be mapped.  This did not affect results.
  es <- con$getGEMatrix(mats, outputType = "raw", annotation = "ImmSig")
  
  # SDY212 has special requirements:
  # 1. Probe with NAs must be taken out before normalization.
  # This is a known issue with ImmPort and is still under investigation.
  # 2. prior to sampleNames mapping, one of duplicated biosamples must be removed.
  if(sdy == "SDY212"){
    # expr matrix
    em <- exprs(es)
    em <- em[ -(grep("ILMN_2137536", rownames(em))), ]
    cnames <- colnames(em)
    rnames <- rownames(em)
    em <- preprocessCore::normalize.quantiles(log2(em))
    colnames(em) <- cnames
    rownames(em) <- rnames
    em <- em[ , -(grep("BS694717.1", colnames(em))) ]
    
    # features
    features <- fData(es)
    features <- features[ -(grep("ILMN_2137536", features$FeatureId)), ]
    rownames(features) <- features$FeatureId
    
    # pheno
    pheno <- pData(es)
    pheno <- pheno[ -(grep("BS694717.1", pheno$biosample_accession)), ]
    
    es <- ExpressionSet(assayData = em,
                        featureData = AnnotatedDataFrame(features),
                        phenoData = AnnotatedDataFrame(pheno))
  }
  
  es <- con$EMNames(es) # map biosample ids to subject ids
  em <- exprs(es)
  
  if(sdy %in% c("SDY63", "SDY404")){
    cnames <- colnames(em)
    rnames <- rownames(em)
    em <- preprocessCore::normalize.quantiles(log2(em)) # log2 MUST be before norm!
    colnames(em) <- cnames
    rownames(em) <- rnames
    
  }else if(sdy == "SDY67"){
    
    # rownames are gene_symbols b/c SDY67 is RNAseq
    # Normalization is done using functions from DESeq
    # according to original analysis.
    countTable <- em
    condition <- colnames(em)
    cds <- newCountDataSet(countTable, condition)
    cds <- estimateSizeFactors(cds) ## estimate size factor
    cdsBlind <- estimateDispersions(cds, method="blind" )
    vsd <- varianceStabilizingTransformation(cdsBlind)
    em <- exprs(vsd)
    
    # gene_set actually has some as lower case, e.g. 'C10orf11', but to
    # reproduce manuscript, must use toupper(), which causes these
    # NOT to be picked up in downstream analysis.
    rownames(em) <- toupper(rownames(em))
  }
  
  # Update sub-ids format: "SUB123123.212_d0" to "SUB123123_d0"
  colnames(em) <- gsub(paste0(".", gsub("SDY", "", sdy)),"",colnames(em), fixed = T)
  em <- data.frame(em, stringsAsFactors = F)
  
  # Not sure why, but the correct "orf" gene symbols are not being found in the fData?!
  # Even though they are present in the schema (microarray > featureAnnotation)
  em$geneSymbol <- if( sdy != "SDY67" ){
    fData(es)[["gene_symbol"]]
  }else{
    rownames(em)
  }
  
  # Remove remaining dup biosample based on Yale HIPCMetaModuleAnalysis.R comments.
  # Remove two subs not found in original file, but needed for normalization.
  if(sdy == "SDY212"){
    # em <- em[ , -(grep("SUB134307.*", colnames(em))) ] #YK
    em <- em[ , -(grep("SUB134242_d0", colnames(em))) ]
    em <- em[ , -(grep("SUB134267_d0", colnames(em))) ]
    
    # Yale collaborators have not yet updated GEO database to reflect these changes.
  }else if( sdy == "SDY404"){
    em <- swapCols(em, "SUB120473_d0", "SUB120474_d0") # Incorrect Sex
    em <- swapCols(em, "SUB120460_d0", "SUB120485_d28") # Strong evidence
  }
  
  # --------GET RAW HAI DATA-----------
  
  # Adjust according to Yuri Kotliarov code.
  # NOTE: SDY80 raw HAI data in ImmuneSpace not available yet, therefore using preloaded.
  # labkey.url.path <- paste0("/Studies/", sdy, "/") # for IS site only
  rawHai <- if(sdy != "SDY80"){ con$getDataset("hai") }else{ SDY80_rawtiterdata_v2 }
  adjHai <- adjust_hai(sdy, rawHai)
  
  # Age needed for figure 2 and SDY80 filtering
  dotNum <- gsub("SDY", ".", sdy, fixed = T)
  demo <- con$getDataset("demographics")
  adjHai$Age <- sapply(adjHai$subject, FUN = function( sub ){
    trg <- demo[ which(demo$participant_id == paste0(sub, dotNum))]
    return( trg$age_reported )
  })
  
  if(sdy == "SDY80"){ adjHai <- adjHai[ which(adjHai$Age < 36), ] }
  
  # Need study id for figure 2
  adjHai$study <- sdy
  
  # ---------PREP + PUSH TO HOLDERS-------
  eset_list[[sdy]] <- prepEset(em, adjHai, d0 = F)
  if( sdy == "SDY80" ){ sdy80all <- prepEset(em, adjHai, d0 = F) }
  
}

# fn.eset = file.path(PROJECT_DIR, "generated_data","HIPC", paste0(sdy, "_IS_eset.rds"))
fn.esets = file.path(PROJECT_DIR, "generated_data","HIPC", "HIPC_IS_esets.rds")
saveRDS(eset_list, fn.esets)

