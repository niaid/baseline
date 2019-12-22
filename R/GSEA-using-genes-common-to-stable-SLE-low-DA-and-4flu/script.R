library(data.table)
library(fgsea)
library(ggplot2)

dn.out = file.path(PROJECT_DIR, "generated_data/GSEA-using-genes-common-to-stable-SLE-low-DA-and-4flu")
dir.create(dn.out, showWarnings = F)
dn.fig = file.path(PROJECT_DIR, "figure_generation/SLE-Sig")
dir.create(dn.fig, showWarnings = F)

WGCNA.mods = fread(file.path(PROJECT_DIR, "generated_data/WGCNA-modules-from-SLE-low-DA/SLE-low-34sbj-9601probes-gene.in.module.minModSize20.signed-hybrid.txt"))
tab0 = fread(file.path(PROJECT_DIR, "generated_data/MetaDE/MetaDE-Effect-Size-from-four-datasets-16083genes.csv"))

common.genes.mask = (tab0$gene %in% WGCNA.mods$Symbol)

tab1 = tab0[common.genes.mask]
stopifnot(all(diff(tab0$zval) <=0))  # already sorted in descending order of zval
stopifnot(all(tab1$gene %in% WGCNA.mods$Symbol))

ranks = tab1$zval
names(ranks) = tab1$gene

mod.names = unique(WGCNA.mods$Module)
mod.names = setdiff(mod.names, "grey")

mods = lapply(mod.names, function(x){WGCNA.mods[which(Module %in% x)]$Symbol})
names(mods) <- mod.names

set.seed(2017072509)
res.fgsea = fgsea(mods, ranks, nperm=100)

fn = file.path(dn.out, "GSEA-enrichment-of-WGCNA-modules-in-4flu-datasets-using-only-genes-in-SLE-dataset.csv")
fwrite(file=fn, res.fgsea)

mod = "brown"
p = plotEnrichment(mods[[mod]], ranks) +
  ggtitle(sprintf("p = %0.3f", res.fgsea$pval[res.fgsea$pathway=="brown"]))

# modify ggplot object
i.param.green = which((sapply(p$layer, function(x)x$aes_params$colour) %>% unlist()) == "green")
for(i in i.param.green) {
  p$layers[[i]]$aes_params$colour = "brown"
}
i.param.red = which((sapply(p$layer, function(x)x$aes_params$colour) %>% unlist()) == "red")
for(i in rev(i.param.red)) {
  p$layers[[i]] <- NULL
}

fn.fig = file.path(dn.fig, sprintf("%s-enrichment.plot-4flu.response.pdf", mod))
ggsave(fn.fig, plot=p, w=6, h=2.5)

