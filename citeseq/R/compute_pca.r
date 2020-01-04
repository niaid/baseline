library(tidyverse)
library(Seurat)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_filtered.rds"
h1 = readRDS(fn)

dn.fig = "figures/TSNE"
dir.create("figures", showWarnings = F)
dir.create(dn.fig, showWarnings = F)

# Run PCA on protein expression
protein = rownames(h1@assay$CITE@raw.data)
protein = protein[-grep("Mouse|Rat",protein)]
h1 <- ScaleData(h1, assay.type = "CITE", display.progress = FALSE, do.par = TRUE, num.cores = 8)
h1 <- RunPCA(h1, assay.type = "CITE", 
             seed.use = 0, pc.genes = protein,
             pcs.print = 0, reduction.name = "protpca")

fn.new = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_PCA.rds"
saveRDS(h1, fn.new)

# elbow plot for PC inclusion in clustering 
y = h1@dr$protpca@sdev
x = 1:20
df = cbind(x,y) %>% base::as.data.frame() 
ggplot(df, mapping = aes(x = x , y = y)) + 
  geom_point() + xlab("PC") + ylab("SD of PC") +
  ggtitle("H1N1 protein day0") + theme_light()

ggsave(file.path(dn.fig, "PCA_elbow.png"))

