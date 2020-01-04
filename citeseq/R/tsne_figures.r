library(Seurat)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

dn.fig = "figures/TSNE"
dir.create(dn.fig, showWarnings = F, recursive = T)

h1 = SetAllIdent(h1, id="batch")

# to check batch effect
p_all = seq(10, 170, by = 40)
for(p in p_all) {
  cat(p, " ")
  DimPlot(h1, reduction.use = paste0("tsne_p",p), pt.size = 0.2) +
    ggtitle(glue::glue("perplexity {p}"))
  ggsave(glue::glue("{dn.fig}/TSNE_perplex.{p}_batches.png"), w=7, h=4)
}

res = c(0.1, 0.3, 1)

# clusters
for(k in 1:(length(res)-1)) {
  cat(k, " ")
  h1 = SetAllIdent(h1, id=paste0("p3_dist_",k))
  for(p in p_all) {
    DimPlot(h1, reduction.use = paste0("tsne_p",p), pt.size = 0.2, cols.use=pals::cols25(length(unique(h1@ident)))) +
      ggtitle(glue::glue("perplexity {p}, res.{res[k]}"))
    ggsave(glue::glue("{dn.fig}/TSNE_perplex.{p}_res.{res[k]}.png"), w=7, h=4)
  }
}


# Plot for Figure
h1 = SetAllIdent(h1, id="K1")
p = 130
# cm = pals::cols25(11)[-6]
cm = pals::glasbey(13)[c(1:3,5:9,13,11)]# %>% as.character()
DimPlot(h1, reduction.use = paste0("tsne_p",p), pt.size = 0.2, 
        cols.use=cm, no.axes = T, no.legend = T)
ggsave(glue::glue("{dn.fig}/TSNE_figure_p{p}.png"), w=7, h=4)
ggsave(glue::glue("{dn.fig}/TSNE_figure_p{p}.pdf"), w=7, h=4)

