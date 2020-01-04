library(corrplot)

dn.fig = "figures/sig_ranked_correlation/"
dir.create(dn.fig, showWarnings = F, recursive = T)

df.scores = fread("results/Score_table.txt")

df.scores.ranked = df.scores %>% 
  mutate_at(vars(matches("\\(")), rank, na.last="keep")

score.mat.ranked = df.scores.ranked %>% 
  dplyr::select(-response) %>% 
  tibble::column_to_rownames("subject") %>% data.matrix()

cc = cor(score.mat.ranked, method = "spearman", use = "pairwise.complete.obs")
ccm = cor.mtest(score.mat.ranked, method = "spearman")
ccp = ccm$p

sname = colnames(score.mat.ranked)
svar = glue::glue("`{sname}`")
xvar = which(sname == "TGSig (microarray)")

for(k in seq_along(sname)[-xvar]) {
  cat(sname[k],"\n")
  ccr = cc[xvar, k]
  ccp = ccm$p[xvar, k]
  df.cc = data.frame(label = glue::glue("r = {format(ccr, digits=2)} p = {format(ccp, digits=2)}"))
  ggplot(df.scores.ranked, aes_string(svar[xvar], svar[k], fill="response")) +
    geom_point(shape = 21, size=3, col="black", show.legend = F) +
    scale_fill_manual(values = c("black","white")) +
    geom_text(data = df.cc, aes(label=label), x = Inf, y=-Inf, hjust=1.1, vjust = -1, col="black", inherit.aes = F) +
    xlab(paste0(sname[xvar], ", ranked")) +
    ylab(paste0(sname[k], ", ranked")) +
    theme_bw()
  fn.fig = glue::glue("{dn.fig}/{sname[k]}_vs_{sname[xvar]}")
  ggsave(paste0(fn.fig, ".png"), w=2.6, h=2.5)
  ggsave(paste0(fn.fig, ".pdf"), w=3, h=2.5, useDingbats=F)
}

