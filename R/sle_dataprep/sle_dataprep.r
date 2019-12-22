library(Biobase)

fn.in = "data/SLE/expression/SLE_Longitudinal_972_eset.RData"
load(file.path(PROJECT_DIR, fn.in), verbose = T)

# eset expression data matrix
dat = exprs(eset)

# eset phenotype data (sample info, demographics, clinical, CBC, SLEDAI components, etc...)
info = pData(eset)
colnames(dat) = info$SAMPLE_NAME

dir.create("./generated_data/SLE", showWarnings = F)

# probe to gene symbol mapping
probes = featureData(eset)@data %>% tibble::rownames_to_column("ID")
probes.map = probes[,c("ID", "ILMN_Gene")]
colnames(probes.map) = c("ID","gene")
probes.map = probes.map[grep("///",probes.map$gene,invert=T),]
probes.map = probes.map[apply(probes.map!="" & !is.na(probes.map), 1, all, na.rm=T),]
fn.map = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_probe_map.txt")
fwrite(probes.map, file=fn.map, sep="\t", quote=F)

# select best probe for each gene
source(file.path(PROJECT_DIR, "R/functions/pick.probeset.r"))
pick.probeset(eset, fn.map)


# log2-transform epxression data
dat = log2(dat)

# output gene expression data and sample info
fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix.txt")
fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info.txt")
dat %>% as.data.frame() %>% tibble::rownames_to_column("ID") %>% 
  fwrite(fn.ge, sep="\t", quote=T)
fwrite(info, fn.si, sep="\t", quote=T)

