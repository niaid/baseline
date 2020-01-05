library(Seurat)
library(tmod)


fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

bg = h1@data %>% rownames()

genes = readRDS("sig/sig.list.RDS")$CD40.act

res = tmodHGtest(genes, bg, mset="LI")
df = res %>% 
  mutate(label = glue::glue("{ID} - {Title}") %>% fct_inorder() %>% fct_rev())

ggplot(df, aes(label, -log10(adj.P.Val))) +
  geom_col() + 
  geom_hline(yintercept = -log10(0.05), col="red", lty=2) +
  xlab("") +
  ylab("adjusted p-value, -log10") +
  coord_flip() +
  theme_classic()

dn.fig = "figures/CD40act"
dir.create(dn.fig, showWarnings = F, recursive = T)

ggsave(file.path(dn.fig, "CD40act_genes_HG.test_p.adj_barplot.png"), w=8, h=4)
ggsave(file.path(dn.fig, "CD40act_genes_HG.test_p.adj_barplot.pdf"), w=8, h=4)
