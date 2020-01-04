library(tidyverse)
library(Seurat)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_PCA.rds"
h1 = readRDS(fn)

fn.new = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE.rds"

for(p in c(seq(10, 170, by = 40))) {
  cat(p, " ")
  h1 = RunTSNE(h1, reduction.use = "protpca", reduction.name = paste0("tsne_p",p), 
               dims.use = 1:7, seed.use = 0, perplexity=p)
  saveRDS(h1, fn.new)
}

