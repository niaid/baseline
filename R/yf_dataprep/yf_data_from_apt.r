fn.ge = file.path(PROJECT_DIR, "data", "YF", "expression", "rma-sketch.summary.txt")
fn.si = file.path(PROJECT_DIR, "data", "YF", "expression", "sample_info.txt")

dat = fread(fn.ge, sep="\t", skip="probeset_id", data.table = F) %>% 
  dplyr::rename(ID=probeset_id)

info = fread(fn.si, sep="\t")
stopifnot(all.equal(colnames(dat)[-1], info$file.name))

info$subject = sub("_Day_\\d+$","",info$sample.name)

info$Response = "middle"
si = info$Trial==1
info$Response[si] = ifelse(info$Neuralizing.Antibody.Titer[si]<160, "low",
                           ifelse(info$Neuralizing.Antibody.Titer[si]>160, "high", "middle"))
si = info$Trial==2
info$Response[si] = ifelse(info$Neuralizing.Antibody.Titer[si]<640, "low",
                           ifelse(info$Neuralizing.Antibody.Titer[si]>640, "high", "middle"))

dir.create("./generated_data/YF", showWarnings = F)
						   
fn.ge.out = file.path(PROJECT_DIR, "generated_data", "YF", "YF_GE_matrix.txt")
fn.si.out = file.path(PROJECT_DIR, "generated_data", "YF", "YF_sample_info.txt")
fwrite(dat, fn.ge.out, sep="\t")
fwrite(info, fn.si.out, sep="\t")


# probe to gene symbol mapping (from hgu133plus2.db package)
# library(hgu133plus2.db)
# probes.map = as.data.frame(hgu133plus2SYMBOL)
# colnames(probes.map) = c("ID","gene")
# fn.map = file.path(PROJECT_DIR, "generated_data", "YF", "YF_probe_map.txt")
# fwrite(probes.map, file=fn.map, sep="\t", quote=F)

# probe to gene symbol mapping (generated from hgu133plus2.db package at the time of analysis)
fn.map = file.path(PROJECT_DIR, "data", "YF/expression/YF_probe_map.txt")
file.copy(from = fn.map, to = file.path(PROJECT_DIR, "generated_data/YF/"))
fn.map = file.path(PROJECT_DIR, "generated_data", "YF/YF_probe_map.txt")

# select best probe for each gene
library(Biobase)
source(file.path(PROJECT_DIR, "R/functions/pick.probeset.r"))
eset = ExpressionSet(dat %>% tibble::column_to_rownames("ID") %>% data.matrix())
pick.probeset(eset, fn.map)
