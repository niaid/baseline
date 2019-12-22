fn.in = file.path(PROJECT_DIR, "generated_data","MetaDE", "MetaDE-Effect-Size-from-four-datasets-16083genes.csv")
meta.res = fread(fn.in)

# load CD38 signature genes
fn.sig = file.path(PROJECT_DIR, "generated_data", "signatures", "CD38_ge_sig.txt")
cd38.genes = fread(fn.sig, header = F) %>% unlist(use.names=F) %>% toupper()

meta.cd38 = meta.res %>% 
  dplyr::filter(toupper(gene) %in% cd38.genes) %>% 
  mutate(mu.ci = qnorm(0.975)*sqrt(mu.var)) %>% 
  mutate(gene = factor(gene, levels=rev(cd38.genes)))

gene.na = (levels(meta.cd38$gene) %in% meta.cd38$gene)+1
gene.clr = c("grey50","black")[gene.na]

ggplot(meta.cd38, aes(gene, mu.hat)) +
  geom_point() +
  geom_errorbar(aes(ymin=mu.hat-mu.ci, ymax=mu.hat+mu.ci), width=0.2) +
  scale_x_discrete(drop=FALSE) +
  xlab("Gene") +
  ylab("Effect size") +
  geom_hline(yintercept=0) +
  coord_flip(ylim = c(-0.45, 1.2)) +
  ggtitle("Influenza vaccination meta-analysis") +
  theme_bw() + theme(axis.text.y = element_text(colour=gene.clr))

dn.fig = file.path(PROJECT_DIR, "figure_generation","MetaDE")
dir.create(dn.fig, showWarnings = F)  
fn.fig = file.path(dn.fig, "MetaDE_flu_forest_plot")
ggsave(paste0(fn.fig, ".png"),w=4, h=3)
ggsave(paste0(fn.fig, ".pdf"),w=4, h=3)

