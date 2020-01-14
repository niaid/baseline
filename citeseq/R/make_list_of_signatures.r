# Make a list of signatures to test

sig.list = list()
sig.list[["TGSig"]] = fread(file.path("../generated_data/signatures/CD38_ge_sig.txt"), header=F) %>% pull(1)
sig.list[["SLE.sig"]] = fread(file.path("../generated_data/fgsea_with_wgcna_modules/brown-leading-edge-isv050-87genes.txt"))$gene
sig.list[["IFN26"]] = fread(file.path("../generated_data/signatures/IFN26_ge_sig.txt"), header=F) %>% pull(1)

# BTMs
library(tmod)
data(tmod)
mod.id.ifn = c("DC.M1.2","DC.M3.4","DC.M5.12")
sig.list[["IFN"]] = tmod$MODULES2GENES[mod.id.ifn] %>% unlist() %>% unique() %>% sort()

mod.id = "LI.M165"
sig.list[[mod.id]] = strsplit(getGenes(mod.id)$Genes, ",")[[1]]

sig.list[["CD40.act"]] = fread("data/B_CD40act_genes_JI2014&Blood2004.txt", header = F) %>% pull(1)

dir.create("sig", showWarnings = F)
saveRDS(sig.list, "sig/sig.list.RDS")


# Make a table

sig.list = readRDS("sig/sig.list.RDS")
sig.list.2 = lapply(sig.list, function(x) paste(x, collapse=", "))
ng = sapply(sig.list, length)
sig.df = stack(sig.list.2) %>% 
  dplyr::select(Signature=ind, Genes=values) %>% 
  add_column(N_genes = ng, .before = "Genes")
fwrite(sig.df, "sig/signature_genes_table.txt", sep="\t")
