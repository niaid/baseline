df.scores = fread("results/Score_table.txt") %>% 
  dplyr::mutate(subject = as.character(subject) %>% factor())
df.scores.rank = df.scores %>% 
  mutate_at(vars(matches("\\(")), rank, na.last="keep")

dn.fig = "figures/sig_ranked_correlation"
dir.create(dn.fig, showWarnings = F, recursive = T)

cc = cor.test(df.scores$`CD38hi (flow)`, df.scores$`CD40act (C3.1.0)`, method = "spearman")
df.txt = data.frame(x=Inf, y=-Inf, label=glue::glue("r = {format(cc$estimate, digits=2)}\np = {format(cc$p.value, digits=2)}"))

ggplot(df.scores.rank, aes(`CD38hi (flow)`, `CD40act (C3.1.0)`, fill=response)) +
  geom_point(shape = 21, size=3, col="black", show.legend = F) +
  scale_fill_manual(values = c("black","white")) +
  geom_text(aes(x, y, label=label), col="black", hjust=1.1, vjust=-0.5, data=df.txt, inherit.aes = F) +
  xlab("CD20+CD38++ B cells (% of parent)") +
  ylab("CD40act score in C3.1.0 cluster") +
  theme_bw()
ggsave(file.path(dn.fig, "CD40act_C3.1.0_vs_CD38hi.flow.png"), w=3.2, h=3)
ggsave(file.path(dn.fig, "CD40act_C3.1.0_vs_CD38hi.flow.pdf"), w=3.2, h=3, useDingbats=F)
