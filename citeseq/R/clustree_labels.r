library(Seurat)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE.rds"
h1 = readRDS(fn)

# labels from clustree (generated manually)
df.label = fread("data/clustree_node_labels_withCellTypeLabels.txt", sep="\t", data.table = F, fill = T)

h1@meta.data = h1@meta.data %>% 
  add_column(K0 = "pseudo-bulk")
  
for(clustering in 1:3) {
  cat("\n", clustering)
  clust.field = glue::glue("p3_dist_{clustering}")
  
  clust.names = h1@meta.data %>% pull(clust.field) #%>% unique() %>% as.numeric() %>% sort()
  df.labels.cl = df.label %>% dplyr::filter(clustering == clust.field)
  h1@meta.data = h1@meta.data %>% 
    add_column(!!(paste0("K",clustering)) := df.labels.cl$label[match(clust.names, df.labels.cl$cluster)])
}

fn.out = str_replace(fn, ".rds", "_labels.rds")
saveRDS(h1, file=fn.out)
