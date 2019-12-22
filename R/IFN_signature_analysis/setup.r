library(tmod)
data(tmod)

# load genes from WGCNA modules from SLE data analysis
fn.wgcna = file.path(PROJECT_DIR, "generated_data", "WGCNA-modules-from-SLE-low-DA", 
                     "SLE-low-34sbj-9601probes-gene.in.module.minModSize20.signed-hybrid.txt")
df.wgcna = fread(fn.wgcna, data.table=F)
genes.wgcna = split(df.wgcna$Symbol, df.wgcna$Module)

mod.id = c("LI.M75", "LI.M150", "LI.M165")
imod = tmod$MODULES$ID %in% mod.id
genes.btm.list = tmod$MODULES2GENES[imod]

genes.btm = unlist(genes.btm.list) %>% unique()

ibrown = genes.wgcna$brown %in% genes.btm
gene.sig = genes.wgcna$brown[ibrown]

dir.create("./generated_data/IFN26", showWarnings = F)
dir.create("./figure_generation/IFN26", showWarnings = F)

fn.sig = file.path("./generated_data/signatures", "IFN26_ge_sig.txt")
gene.sig %>% as.data.frame() %>% 
  fwrite(fn.sig, col.names = F)
