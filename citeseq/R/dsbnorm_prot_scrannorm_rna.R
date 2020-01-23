suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

source("R/functions/dsb_normalization_functions.R")
dir.create("data/normalization_data")


# read in final singlets make list of seurat objects indexed by batch. 
h1 = readRDS(file = "data/H1_day0_demultilexed_singlets.RDS") %>% SetAllIdent(id = "batch")
h1b1 = SubsetData(h1, ident.use = "1")
h1b2 = SubsetData(h1, ident.use = "2")

## Protein normalization with DSB normalization
# load object with all cells and subset the negative cells,
s = readRDS(file = "data/neg_control_object.rds")
batch_vector = unique(s@meta.data$batch)
s = SetAllIdent(s, id = "batch")
neg_adt = list()
for (i in 1:length(batch_vector)) {
  neg_adt[[i]] = s %>% SetAllIdent(id = "batch") %>% SubsetData(ident.use = batch_vector[i])
}
rm(s) 
gc()

# filter by nGene 
# get list of protein matrix by batch for negative control drops.  
neg_adt = lapply(neg_adt, function(x){ SubsetData(x, accept.high = 80, subset.name = "nGene") })
neg_adt = lapply(neg_adt, function(x){ x@assay$CITE@raw.data })

# make list of positive protein matrices by batch 
stained = list(h1b1, h1b2)
pos_adt = lapply(stained, function(x){x@assay$CITE@raw.data})

# apply denoised scaled by background protein normalization. 
dsb_norm = list()
for (i in 1:length(neg_adt)) {
  dsb_norm[[i]] = 
    DSBNormalizeProtein(cell.columns.protein.matrix = pos_adt[[i]], control.protein.matrix = neg_adt[[i]],
                        define.pseudocount = TRUE, pseudocount.use = 1, denoise_counts = TRUE, 
                        isotype.control.name.vec =  c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", 
                                                      "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT"))
  
}

# these data get added to object by batch at end of script
saveRDS(dsb_norm, file = "data/normalization_data/dsb_norm_list_bybatch.rds")


##mRNA normalization with scran 
# normalize by batch, use multibatch norm on h1 . 
sc = list(h1b1, h1b2)
sc =  lapply(sc, function(x){ Convert(from = x, to = "sce") })  
suppressMessages(library(scater))
sc = lapply(sc, FUN = calculateQCMetrics)


# outlier cells based on lib size n median absolute deviations > or < median lib size 
low = lapply(sc, function(x){ isOutlier(x$log10_total_counts, type = "lower", nmads = 3) })
high = lapply(sc, function(x){ isOutlier(x$log10_total_counts, type = "higher", nmads = 3) })

# remove outlier cells from the matrix. 
outlier_cell = list()
for (i in 1:length(sc)) {
  outlier_cell[[i]] = low[[i]] | high[[i]]
  sc[[i]] = sc[[i]][ ,!outlier_cell[[i]]]
}

#normalize with multibatch normalization (now part of the batchelor package ) 
library(scran)
sc = lapply(sc, function(x){ computeSumFactors(x, min.mean = 0.1 , get.spikes = FALSE) })
batch.norm = multiBatchNorm(sc[[1]],  sc[[2]], min.mean = 0.1)
saveRDS(batch.norm, file = "data/normalization_data/sce_list_multibatchnormalized.rds")


## Add data back to Seurat object, subset by retained non outlier cells. 
b1norm = logcounts(batch.norm[[1]])
b2norm = logcounts(batch.norm[[2]])

# normalized mRNA 
h1_norm = cbind(b1norm, b2norm)
vars_add = c("total_features_by_counts", "log10_total_features_by_counts", "total_counts",
             "log10_total_counts", "pct_counts_in_top_50_features")
# metadata 
b1meta = 
  colData(batch.norm[[1]]) %>% 
  as.data.frame() %>% 
  dplyr::select(vars_add)
b1meta$sizefactors = sizeFactors(batch.norm[[1]])
b2meta = 
  colData(batch.norm[[2]]) %>% 
  as.data.frame() %>% 
  dplyr::select(vars_add)
b2meta$sizefactors = sizeFactors(batch.norm[[2]])
h1_meta_add = rbind(b1meta, b2meta)
## Add back to seurat object 
h1 = h1 %>% 
  SubsetData(cells.use = rownames(h1_meta_add)) %>% 
  AddMetaData(metadata = h1_meta_add) %>% 
  SetAssayData(new.data = h1_norm, assay.type = "RNA", slot = "data")
# tell seurat the data is normalized w a custom method. 
h1@calc.params$NormalizeData = list(normalization.method = "scran")


## Finally, Add the DSB normalized CITE assay to the object with scrannormalized RNA counts. 
# dsb_norm = readRDS("data/dsb_norm_list_bybatch.rds") 
# subset mRNA outlier cells from adt data 
for (i in 1:length(dsb_norm)) {
  dsb_norm[[i]] = dsb_norm[[i]][ ,batch.norm[[i]]$barcode_check]
}
h1_CITE = base::do.call(cbind, dsb_norm[c(1,2)])
h1 = SetAssayData(h1, new.data = h1_CITE, assay.type = "CITE", slot = "data")

# save the processed data with scran normalized mRNA and dsb normalized protein 
saveRDS(h1, file = "data/H1_day0_scranNorm_adtbatchNorm.rds")

