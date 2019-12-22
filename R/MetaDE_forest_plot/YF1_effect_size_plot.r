fn.in = file.path(PROJECT_DIR, "generated_data","MetaDE", "YF1-hedges-g.csv")
yf.res = fread(fn.in)

# load CD38 signature genes
fn.sig = file.path(PROJECT_DIR, "generated_data", "signatures", "CD38_ge_sig.txt")
cd38.genes = fread(fn.sig, header = F) %>% unlist(use.names=F) %>% toupper()

yf.cd38 = yf.res %>% 
  dplyr::filter(toupper(Gene) %in% cd38.genes) %>%
  mutate(Gene = factor(Gene, levels=rev(cd38.genes)))

gene.na = (levels(yf.cd38$Gene) %in% yf.cd38$Gene)+1
gene.clr = c("grey50","black")[gene.na]

ggplot(yf.cd38, aes(Gene, Hedges.g)) +
  geom_point() +
  geom_errorbar(aes(ymin=Conf95.lower, ymax=Conf95.upper), width=0.2) +
  scale_x_discrete(drop=FALSE) +
  ylab("Effect size") +
  geom_hline(yintercept=0) +
  coord_flip(ylim = c(-1.6, 2.3)) + 
  ggtitle("Yellow Fever vaccination (trial #1)") +
  theme_bw() + theme(axis.text.y = element_text(colour=gene.clr))

dn.fig = file.path(PROJECT_DIR, "figure_generation","MetaDE")
dir.create(dn.fig, showWarnings = F)  
fn.fig = file.path(dn.fig, "YF_effect_size_plot")
ggsave(paste0(fn.fig, ".png"),w=4, h=3)
ggsave(paste0(fn.fig, ".pdf"),w=4, h=3)

