# generate gene list from genes highly correlated with CD38-high B cells  

source("R/functions/load_sig.r")

dir.create(file.path(PROJECT_DIR, "generated_data", "signatures"), showWarnings = F)

fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_genes.txt")
gene.sig = load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=10)
fn.cd38.sig = file.path(PROJECT_DIR, "generated_data", "signatures", "CD38_ge_sig.txt")
gene.sig %>% as.data.frame() %>% 
  fwrite(fn.cd38.sig, col.names = F)


# generate gene lists from BTMs related to plasma cells  

library(tmod)
data(tmod)
mod.id.li = tmod$MODULES %>% dplyr::filter(grepl("plasma",Title, ignore.case = F), Category=="immune") %>% 
  dplyr::select(ID) %>% 
  unlist(use.names=F)
mod.id.dc = c("DC.M4.11","DC.M7.7","DC.M7.32")
mod.id = c(mod.id.dc, mod.id.li)
mod.gene = tmod$MODULES2GENES[mod.id]

for(k in mod.id) {
  fn.pb.sig = file.path(PROJECT_DIR, "generated_data", "signatures", 
                        sprintf("PB_%s_ge_sig.txt", k))
  mod.gene[[k]] %>% as.data.frame() %>% 
    fwrite(fn.pb.sig, col.names = F)
}

