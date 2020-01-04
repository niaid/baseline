fn.info = file.path("../generated_data/CHI/CHI_sample_info_2_CD38hi.txt")
df.info = fread(fn.info) %>% 
  dplyr::filter(time==0, Response %in% c("low","high"))

sigs = c("IFN26", "LI.M165")
sig.names = c("IFN-I-DCact", "LI.M165 (DC activation module)") %>% setNames(sigs)
clustering = 1
clust.name = "C9"

dn.fig = "figures/sig_test_boxplots"
dir.create(dn.fig, showWarnings = F, recursive = T)

for(s in sigs) {
  df.mat = fread(glue::glue("results/sig_scores/scores_{s}_{clustering}.txt"))
  
  # add Responders
  df = df.mat %>% left_join(df.info %>% dplyr::select(subject, Response), by="subject") %>% 
    mutate(Response = factor(Response, levels=c("low","high"))) %>% 
    dplyr::select(subject, Response, C9)
  
  X = df %>% pull(clust.name)
  Y = df$Response
  tt = t.test(X[Y=="high"], X[Y=="low"], data=df, alternative="greater")
  wt = wilcox.test(X[Y=="high"], X[Y=="low"], data=df, alternative="greater")
  
  df.txt = data.frame(x=-Inf, y=Inf, label=glue::glue("p = {format(wt$p.value, digits=2)}"))
  
  cm = pals::glasbey(13)[c(1:3,5:9,13,11)]
  yrng = c(-1,1.55)
  dotsz = 2
  ggplot(df, aes(Response, C9, group=Response)) +
    coord_flip(ylim = c(-0.55, 0.65)) + 
    geom_boxplot(width=0.6, fill=cm[10], position = position_dodge2(preserve = "total"), lwd = 0.3, alpha = 0.5, outlier.shape = NA, show.legend = F) +
    geom_dotplot(aes(fill = Response), binaxis = "y", stackdir = "center", dotsize = dotsz, binwidth=0.03,
                 position=position_dodge(), show.legend = F) +
    geom_text(aes(x, y, label=label), col="red", hjust=1, vjust=-0.5, data=df.txt, inherit.aes = F) +
    scale_fill_manual(values = c("white", "black")) +
    ylab(glue::glue("{sig.names[s]} score")) +
    theme_classic()

  fn.fig = glue::glue("{dn.fig}/{s}_C9_score_boxplot")
  ggsave(paste0(fn.fig, ".png"), w=2.8, h=1.5)
  ggsave(paste0(fn.fig, ".pdf"), w=2.8, h=1.5, useDingbats=F)
}
