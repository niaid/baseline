library(fgsea)

source(file.path(PROJECT_DIR, "R/functions/load_sig.r"))

dn.out = file.path(PROJECT_DIR, "generated_data", "fgsea_with_wgcna_modules")
dir.create(dn.out, showWarnings = F)
dn.fig = file.path(PROJECT_DIR, "figure_generation", "SLE-Sig")
dir.create(dn.fig, showWarnings = F)

# load CD38 signature genes
fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_all.genes.txt")
df.cd38.cor = fread(fn.cd38.cor, data.table=F)

# load gene stability
fn.stab = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_genes_stability.txt")
ge.stab = fread(fn.stab, data.table=F)

# load genes from WGCNA modules from SLE data analysis
fn.wgcna = file.path(PROJECT_DIR, "generated_data", "WGCNA-modules-from-SLE-low-DA", 
                     "SLE-low-34sbj-9601probes-gene.in.module.minModSize20.signed-hybrid.txt")
df.wgcna = fread(fn.wgcna, data.table=F)
genes.wgcna = split(df.wgcna$Symbol, df.wgcna$Module)


df.comb = left_join(df.cd38.cor, ge.stab, by="gene") %>% 
  arrange(-cor.mean.sd.ratio)


# GSEA for genes ranked yb correlation with CD20+CD38++ B cell freq. at different ISV cutoff 

isv.seq = seq(0, 0.75, by=0.25)
isv.rank = df.comb$ISV
names(isv.rank) = df.comb$gene


l = list()
set.seed(123)
cat("ISV cutoff: ")
for(i in seq_along(isv.seq)) {
  cat(isv.seq[i],"")
  gi = df.comb$ISV >= isv.seq[i]
  cc.rank = df.comb$cor.mean.sd.ratio[gi]
  names(cc.rank) = df.comb$gene[gi]
  l[[paste0("ISV.",isv.seq[i])]] <- fgsea(genes.wgcna, cc.rank, nperm=10000, maxSize=500)
}

fn.res = file.path(dn.out, "WGCNA_genes_fgsea_results.rds")
saveRDS(l, fn.res)

# fn.res = file.path(dn.out, "WGCNA_genes_fgsea_results.rds")
# l = readRDS(fn.res)


# Enrichment plot for the Brown module at ISV cutoff 0.5 ---

mset = "brown"
isv.th = 0.5
gi = df.comb$ISV >= isv.th
cc.rank = df.comb$cor.mean.sd.ratio[gi]
names(cc.rank) = df.comb$gene[gi]

p = plotEnrichment(genes.wgcna[[mset]], cc.rank) +
  ggtitle(sprintf("p = %0.4f", l$ISV.0.5$pval[l$ISV.0.5$pathway=="brown"]))

# modify ggplot object
i.param.green = which((sapply(p$layer, function(x)x$aes_params$colour) %>% unlist()) == "green")
for(i in i.param.green) {
  p$layers[[i]]$aes_params$colour = "brown"
}
i.param.red = which((sapply(p$layer, function(x)x$aes_params$colour) %>% unlist()) == "red")
for(i in rev(i.param.red)) {
  p$layers[[i]] <- NULL
}
fn.fig = glue::glue("{dn.fig}/{mset}_enrichment.plot_CD38.cor.ranked_ISV.{isv.th}.pdf")
ggsave(fn.fig, plot=p, w=6, h=2.5)


# get leading genes for a module ------------------------------------------
gs = "brown"
le = list()
df.le = data.frame()
for(i in 1:length(l)) {
  le[i] = l[[i]] %>% 
    dplyr::filter(pathway==gs) %>% 
    pull(leadingEdge)
  tmp = data.frame(gene = le[[i]], isv=names(l)[i])
  df.le = bind_rows(df.le, tmp)
}

fn.out = file.path(dn.out, sprintf("CD38_ge_cor_%s.module_leadingEdge_genes.txt", gs))
fwrite(df.le, fn.out, sep="\t", quote = T)

fn.out = file.path(dn.out, sprintf("brown-leading-edge-isv050-87genes.txt", gs))
fwrite(df.le %>% dplyr::filter(isv %in% "ISV.0.5") %>% dplyr::select(gene), fn.out, sep="\t", quote = T)

