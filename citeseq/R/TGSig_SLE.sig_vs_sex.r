library(Seurat)
source("R/functions/gg_color_hue.r")

dn.fig = "figures/male_vs_female"
dir.create(dn.fig, showWarnings = F, recursive = T)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

df.subj = h1@meta.data %>% dplyr::mutate(response = str_remove(adjmfc.time, "d0 ")) %>% 
  dplyr::select(subject=sampleid, response) %>% 
  distinct() %>% 
  dplyr::arrange(subject)

fn.info = "../generated_data/CHI/CHI_sample_info_2_CD38hi.txt"
info = fread(fn.info) %>% dplyr::mutate(subject = as.character(subject)) %>% 
  dplyr::filter(time==0)

fn.tgsig = "results/sig_scores/scores_TGSig_0.txt"
df.tgsig = fread(fn.tgsig) %>% 
  dplyr::rename(TGSig = `pseudo-bulk`) %>% 
  dplyr::mutate(subject = as.character(subject))

fn.sle = "results/sig_scores/scores_SLE.sig_0.txt"
df.sle = fread(fn.sle) %>% 
  dplyr::rename(SLE.sig = `pseudo-bulk`) %>% 
  dplyr::mutate(subject = as.character(subject))

df = df.subj %>% 
  left_join(info, by="subject") %>% 
  left_join(df.tgsig, by="subject") %>% 
  left_join(df.sle, by="subject") %>% 
  dplyr::mutate(gender = factor(gender)) %>% 
  dplyr::mutate(Response = factor(Response, levels=c("low","high"))) %>% 
  gather("Sig","Score", c(TGSig, SLE.sig)) %>% 
  dplyr::mutate(Sig = fct_inorder(Sig))


df.w = df %>% 
  group_by(Sig) %>% 
  do(broom::tidy(wilcox.test(Score ~ gender, data=., exact=F, 
                             paired=F, alternative = "two"))) %>% 
  ungroup() %>% 
  dplyr::mutate(label = glue::glue("p = {format(p.value, digits=2)}"))

clr = gg_color_hue(2)

ggplot(df, aes(gender, Score)) +
  geom_boxplot(aes(fill=gender), alpha=1, outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
  # geom_jitter(aes(col=Response), width = 0.2, height = 0, size=2) +
  facet_wrap(~Sig, nrow=1) +
  scale_fill_manual(values=c(clr[1], "black","white", clr[2])) +
  xlab("Gender") +
  geom_text(data=df.w, aes(label=label), x=1.5, y=Inf, hjust=0.5, vjust = 1.1, size=3, inherit.aes = F) +
  theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
        axis.ticks = element_blank())
ggsave(file.path(dn.fig, "TGSig_SLE.sig_vs_sex.png"), w=4, h=4)
ggsave(file.path(dn.fig, "TGSig_SLE.sig_vs_sex.pdf"), w=4, h=4, useDingbats=F)
