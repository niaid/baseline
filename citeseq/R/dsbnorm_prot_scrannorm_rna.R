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
  select(vars_add)
b1meta$sizefactors = sizeFactors(batch.norm[[1]])
b2meta = 
  colData(batch.norm[[2]]) %>% 
  as.data.frame() %>% 
  select(vars_add)
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


sessionInfo()
# MPM -- ran on http://ai-rstudioprd1.niaid.nih.gov:8787/
# R version 3.5.2 (2018-12-20)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux Server 7.6 (Maipo)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] scran_1.10.2                scater_1.9.20               SingleCellExperiment_1.3.10
# [4] SummarizedExperiment_1.11.6 DelayedArray_0.7.37         BiocParallel_1.15.11       
# [7] matrixStats_0.54.0          Biobase_2.41.2              GenomicRanges_1.33.13      
# [10] GenomeInfoDb_1.17.1         IRanges_2.15.17             S4Vectors_0.19.19          
# [13] BiocGenerics_0.27.1         mclust_5.4.3                magrittr_1.5               
# [16] bindrcpp_0.2.2              here_0.1                    Seurat_2.3.4               
# [19] Matrix_1.2-15               cowplot_0.9.4               forcats_0.3.0              
# [22] stringr_1.4.0               dplyr_0.7.8                 purrr_0.3.2                
# [25] readr_1.3.1                 tidyr_0.8.3                 tibble_2.1.1               
# [28] ggplot2_3.1.1               tidyverse_1.2.1            
# 
# loaded via a namespace (and not attached):
#   [1] reticulate_1.10          R.utils_2.7.0            tidyselect_0.2.5        
# [4] htmlwidgets_1.3          grid_3.5.2               trimcluster_0.1-2.1     
# [7] Rtsne_0.15               munsell_0.5.0            codetools_0.2-15        
# [10] ica_1.0-2                statmod_1.4.30           withr_2.1.2             
# [13] colorspace_1.4-1         knitr_1.22               rstudioapi_0.9.0        
# [16] ROCR_1.0-7               robustbase_0.93-3        dtw_1.20-1              
# [19] gbRd_0.4-11              Rdpack_0.9-0             labeling_0.3            
# [22] lars_1.2                 GenomeInfoDbData_1.1.0   bit64_0.9-7             
# [25] rhdf5_2.25.9             rprojroot_1.3-2          xfun_0.6                
# [28] diptest_0.75-7           R6_2.4.0                 ggbeeswarm_0.6.0        
# [31] locfit_1.5-9.1           hdf5r_1.0.0              flexmix_2.3-14          
# [34] bitops_1.0-6             assertthat_0.2.1         SDMTools_1.1-221        
# [37] scales_1.0.0             nnet_7.3-12              beeswarm_0.2.3          
# [40] gtable_0.3.0             npsurv_0.4-0             rlang_0.3.1             
# [43] splines_3.5.2            lazyeval_0.2.2           acepack_1.4.1           
# [46] broom_0.5.0              checkmate_1.9.1          yaml_2.2.0              
# [49] reshape2_1.4.3           modelr_0.1.2             backports_1.1.3         
# [52] Hmisc_4.2-0              tools_3.5.2              gplots_3.0.1.1          
# [55] RColorBrewer_1.1-2       proxy_0.4-22             dynamicTreeCut_1.63-1   
# [58] ggridges_0.5.1           Rcpp_1.0.0               plyr_1.8.4              
# [61] base64enc_0.1-3          zlibbioc_1.27.0          RCurl_1.95-4.11         
# [64] rpart_4.1-13             pbapply_1.3-4            viridis_0.5.1           
# [67] zoo_1.8-3                haven_1.1.2              cluster_2.0.7-1         
# [70] data.table_1.12.0        lmtest_0.9-36            RANN_2.6.1              
# [73] mvtnorm_1.0-8            fitdistrplus_1.0-14      hms_0.4.2               
# [76] lsei_1.2-0               evaluate_0.13            readxl_1.1.0            
# [79] gridExtra_2.3            compiler_3.5.2           KernSmooth_2.23-15      
# [82] crayon_1.3.4             R.oo_1.22.0              htmltools_0.3.6         
# [85] segmented_0.5-3.0        Formula_1.2-3            snow_0.4-2              
# [88] lubridate_1.7.4          MASS_7.3-51.1            fpc_2.1-11.1            
# [91] cli_1.1.0                R.methodsS3_1.7.1        gdata_2.18.0            
# [94] metap_1.0                bindr_0.1.1              igraph_1.2.4            
# [97] pkgconfig_2.0.2          foreign_0.8-71           xml2_1.2.0              
# [100] foreach_1.4.4            vipor_0.4.5              XVector_0.21.3          
# [103] bibtex_0.4.2             rvest_0.3.2              digest_0.6.18           
# [106] tsne_0.1-3               rmarkdown_1.10           cellranger_1.1.0        
# [109] htmlTable_1.13.1         edgeR_3.23.3             DelayedMatrixStats_1.3.8
# [112] kernlab_0.9-27           gtools_3.8.1             modeltools_0.2-22       
# [115] nlme_3.1-137             jsonlite_1.6             Rhdf5lib_1.3.3          
# [118] BiocNeighbors_1.0.0      viridisLite_0.3.0        limma_3.37.4            
# [121] pillar_1.3.1             lattice_0.20-38          httr_1.3.1              
# [124] DEoptimR_1.0-8           survival_2.43-3          glue_1.3.1              
# [127] png_0.1-7                prabclus_2.2-6           iterators_1.0.10        
# [130] bit_1.1-14               class_7.3-14             stringi_1.4.3           
# [133] HDF5Array_1.9.15         mixtools_1.1.0           doSNOW_1.0.16           
# [136] latticeExtra_0.6-28      caTools_1.17.1.1         irlba_2.3.3             
# [139] ape_5.2   
