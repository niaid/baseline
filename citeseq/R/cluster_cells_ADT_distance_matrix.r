suppressMessages(library(Seurat)) 

fn = "data/H1_day0_scranNorm_adtbatchNorm.rds"
h1 = readRDS(fn)


# specify levels of protein to cluster, remove proteins that didnt stain. 
prot = rownames(h1@assay$CITE@data)
prot_subset = prot[-c(83:86)]


# Subset normalized protein 3 levels. 
h1_adt = GetAssayData(h1, assay.type = "CITE", slot = "data")


# subset used markers 
h1_adt3 = h1_adt[prot_subset, ]


# get distance matrix 
p3_dist = dist(t(h1_adt3))
p3_dist = as.matrix(p3_dist)


# cluster 
res = c(0.1, 0.3, 1, 3.0)
for (i in 1:length(res)) {
  cat(res[i]," ")
  h1 = FindClusters(h1, 
  					distance.matrix = p3_dist,
  					k.param = 50,
                    print.output = F, 
                    resolution = res[i], 
                    random.seed = 1,
                    algorithm = 3,
                    modularity.fxn = 1)
  h1 = StashIdent(h1, save.name = paste0("p3_dist_",i))
  
  saveRDS(h1@meta.data, "data/H1N1_full_clustered_p3dist_multires_metadata.rds")
}

# Add clusters to the original object. (The new object stores the distance matrix and is very large to be saved.)

h1_md = h1@meta.data
h1 = readRDS(fn) # reread the original object
h1 = AddMetaData(h1, metadata = h1_md)
fn.new = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered.rds"
saveRDS(h1, fn.new)

