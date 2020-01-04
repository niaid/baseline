# Calculate average normalized protein expression by cluster**

suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

dn.fig = file.path("figures/cluster_annotation")
dir.create(dn.fig, showWarnings = F, recursive = T)
dir.create("results", showWarnings = F)

fn = file.path("data/H1_day0_scranNorm_adtbatchNorm_dist_clustered.rds")
h1 = readRDS(fn)
md = h1@meta.data %>% dplyr::select(starts_with("p3_dist")) %>% mutate_all(.funs = as.numeric)

# generate cluster table for annotation
df.clust = md[,1:3] %>% 
  gather("clustering","cluster") %>% 
  distinct() %>% 
  arrange(clustering, cluster)
# use this table to manually add cluster label and cell type annotation
fwrite(df.clust, file.path("results/cluster_nodes.txt"), sep="\t")

# relation between clusters from different resolutions
library(clustree)
clustree(md, prefix = "p3_dist_")
# use this figure to label the clusters
ggsave(file.path(dn.fig, "clustree_first.png"), w=12, h=7)


# For each resolution used in clustering, calculate the average expression. 

adt = h1@assay$CITE@data %>% 
  t %>% 
  as.data.frame() %>% 
  rownames_to_column("cell")
proteins = names(adt)[-1]

adt = cbind(adt, md)

# res 0.1 
adt1 = adt %>% 
  group_by(cluster=p3_dist_1) %>% 
  summarize_at(.vars = proteins, .funs = base::mean) %>% 
  tibble::column_to_rownames("cluster") %>% data.matrix()
colnames(adt1) = str_sub(colnames(adt1), start = 1, end = -6)

# res 0.3 
adt2 = adt %>% 
  group_by(cluster=p3_dist_2) %>% 
  summarize_at(.vars = proteins, .funs = base::mean) %>% 
  tibble::column_to_rownames("cluster") %>% data.matrix()
colnames(adt2) = str_sub(colnames(adt2), start = 1, end = -6)

# Res 1.0 
adt3 = adt %>% 
  group_by(cluster=p3_dist_3) %>% 
  summarize_at(.vars = proteins, .funs = base::mean) %>% 
  tibble::column_to_rownames("cluster") %>% data.matrix()
colnames(adt3) = str_sub(colnames(adt3), start = 1, end = -6)

# Res 3.0 
adt4 = adt %>% 
  group_by(cluster=p3_dist_4) %>% 
  summarize_at(.vars = proteins, .funs = base::mean) %>% 
  tibble::column_to_rownames("cluster") %>% data.matrix()
colnames(adt4) = str_sub(colnames(adt4), start = 1, end = -6)


adt.mat = rbind(adt1, adt2, adt3, adt4)


# proteins to use
prot_use = c("CD2","CD3","CD5","CD7","CD4","CD8","CD62L","CD45RA","CD45RO","CD27","CD28","CD278 ", "CD25", 
             "CD127","CD161", "KLRG1","CD195","CD314 ", "CD194","CD103", "CD56","CD57","CD244",
             "CD16","CD14", "CD11b","CD11c", "CD1d", "CD33", "CD13", "CD31", "CD64", "CD163","CD86", "HLA-DR", 
             "CD123", "CD141", "CD71", "CD303", "CD117","CD34", "CD38","CD39", "CD1c", "CD32", 
             "IgM","IgD","IgA","CD19","CD20","CD21","CD24","CD40","CD185","CD196","BTLA")

i.use = match(prot_use, colnames(adt.mat))


# Heatmap

library(ComplexHeatmap)
library(circlize)

clust.level = c(rep(1, nrow(adt1)), rep(2, nrow(adt2)), rep(3, nrow(adt3)), rep(4, nrow(adt4)))
cu = RColorBrewer::brewer.pal(8, name = "OrRd") 
hm = Heatmap(adt.mat[,i.use], name="mean", cluster_rows = T, cluster_columns = F,
             col=cu,
             split =clust.level, cluster_row_slices = F,
             row_names_max_width = max_text_width( rownames(adt.mat), gp = gpar(fontsize = 12) ))
fn.fig = file.path(glue::glue("{dn.fig}/adt_mean_clusters_4levels_Heatmap"))
png(paste0(fn.fig, ".png"), w=1100, h=1100)
draw(hm)
dev.off()

