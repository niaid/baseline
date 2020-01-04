library(Seurat)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered.rds"
h1 = readRDS(fn)

# Remove the clusters that we found to contain mostly doublet cells

i.rm = h1@meta.data$p3_dist_4 %in% c("34","36","38","39")
sum(i.rm)
h1f = SubsetData(h1, cells.use = h1@cell.names[!i.rm])

fn.new = sub(".rds", "_filtered.rds", fn)
saveRDS(h1f, fn.new)

