library(fgsea)
library(tmod)
source(file.path(PROJECT_DIR, "R/functions/load_sig.r"))

# load CD38 signature genes
fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_all.genes.txt")
df.cd38.cor = fread(fn.cd38.cor, data.table=F)

# load gene stability
fn.stab = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_genes_stability.txt")
ge.stab = fread(fn.stab, data.table=F)

df.comb = left_join(df.cd38.cor, ge.stab, by="gene") %>% 
  arrange(-cor.mean.sd.ratio)

# load genes from BTM modules
data(tmod)

for (mset in c("LI","DC")) {
  cat("\nProcessing", mset, "BTM gene sets. ISV cutoff: ")
imod = tmod$MODULES$SourceID==mset
genes.btm = tmod$MODULES2GENES[imod]


isv.seq = seq(0, 0.75, by=0.25)

isv.rank = df.comb$ISV
names(isv.rank) = df.comb$gene

# load("SLE/SLE.RData", verbose = T)
# sle.genes = rownames(SLE$dat)
# isle = toupper(df.comb$gene) %in% toupper(sle.genes)
# sum(isle)

l = list()
set.seed(123)
for(i in seq_along(isv.seq)) {
  cat(isv.seq[i],"")
  gi = df.comb$ISV >= isv.seq[i]
  cc.rank = df.comb$cor.mean.sd.ratio[gi]
  names(cc.rank) = df.comb$gene[gi]
  l[[paste0("ISV.",isv.seq[i])]] <- fgsea(genes.btm, cc.rank, nperm=10000, maxSize=500)
}

dn = file.path(PROJECT_DIR, "generated_data", "fgsea_with_btm_modules")
dir.create(dn, showWarnings = F)
fn.res = file.path(dn, sprintf("BTM.%s_genes_fgsea_results.rds", mset))
saveRDS(l, fn.res)
}

