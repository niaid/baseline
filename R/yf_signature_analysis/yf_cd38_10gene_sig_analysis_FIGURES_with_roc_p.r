library(pROC)
library(cowplot)
source(file.path(PROJECT_DIR, "R/functions/gg_color_hue.r"))

fn.si = file.path(PROJECT_DIR, "generated_data", "YF", "YF_cd38_ge_sig_score.txt")
info = fread(fn.si)%>% 
  dplyr::filter(Response %in% c("low","high")) %>% 
  mutate(Response = factor(Response, levels=c("low","high")),
         label = ifelse(Trial==1, "Yellow Fever", paste0("Yellow Fever, Trial ", Trial)))

roc.df = data.frame()
df.text = data.frame()
N_perm = 1000

for (trial in 1:2) {
  tinfo = info %>% dplyr::filter(Trial == trial)
  
  r = roc(Response ~ CD38_score, data=tinfo, direction="<", quiet = T)
  set.seed(123)
  r.null = vector("numeric", N_perm)
  for(i in 1:N_perm) {
    r.null[i] = roc(Response~sample(CD38_score), data=tinfo, direction="<", quiet = T)$auc
  }
  r.p = sum(r.null > r$auc) / N_perm
  
  r.df = data.frame(Specificity = r$specificities, Sensitivity=r$sensitivities, Trial=trial) %>%
    arrange(Sensitivity)
  
  roc.df = rbind(roc.df, r.df)
  
  
  X = tinfo$CD38_score
  Y = tinfo$Response
  w.pv = wilcox.test(X[Y=="high"], X[Y=="low"], exact=F, alternative="greater")$p.value

  df.text = rbind(df.text, 
                  data.frame(Trial=trial, label.auc=sprintf("AUC = %.2f\np = %.2g", r$auc, r.p), 
                             label.w.pv=sprintf("p = %.2g",w.pv)))
}

roc.df = roc.df %>% 
  mutate(Trial=factor(Trial),
         label = paste0("Yellow Fever, Trial ", Trial)) %>% 
  arrange(Sensitivity, Specificity)
df.text$x = 0.25
df.text$y = 0.1

# df.text$x = c(0.5, 0.5)
# df.text$y = c(roc.df$Sensitivity[which(roc.df$Trial==1 & roc.df$Specificity<=0.5)[1]],
#               roc.df$Sensitivity[which(roc.df$Trial==2 & roc.df$Specificity<=0.5)[1]])
df.text = df.text %>% 
  mutate(label = ifelse(Trial==1, "Yellow Fever", paste0("Yellow Fever, Trial ", Trial)))

clr = gg_color_hue(4) %>% rev()

for (t in 1:2) {
  df.text.2 = df.text %>% 
    dplyr::filter(Trial==t)
  
  p1 = ggplot(info %>% dplyr::filter(Trial==t), aes(x=Response, y=CD38_score, group=Response)) +
    geom_boxplot(aes(fill=label), alpha=0.5, outlier.colour = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
    facet_wrap(~label, nrow=1) +
    scale_fill_manual(values=c("black", "white", clr[4]), name="Response",
                      breaks=c(0,2), labels=c("low","high")) +
    geom_text(data = df.text.2, aes(label=label.w.pv), x=1.5, y=Inf,vjust=1.1, hjust=0.5, size=4, inherit.aes = F) +
    xlab("Response") + ylab("Baseline signature score") +
    coord_cartesian(xlim=c(0.8, 2.2)) +
    theme_bw() + theme(legend.position="none") +
    theme(panel.border = element_blank(), strip.background = element_blank(), 
          strip.text.x = element_text(size=12),
          panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
          axis.ticks = element_blank())
  
  p2 = ggplot(roc.df %>% dplyr::filter(Trial==t), aes(x=Specificity,y=Sensitivity)) + geom_line(size=1) +
    geom_abline(slope=1, intercept=1, lty=2, col="black") +
    scale_x_reverse() + 
    geom_text(data=df.text.2, aes(x=x,y=y,label=label.auc), vjust=-0.5) +
    scale_color_manual(values=c("red","blue")) +
    coord_fixed() + theme_bw() + guides(color="none")
  
  # bottom_row <- do.call("plot_grid", c(p2, align = "h", nrow=1))
  plot_grid(p1, p2, ncol=1, rel_heights = c(2, 1))
  
  fn.fig = file.path(PROJECT_DIR, "figure_generation", sprintf("YF_Trial%d_CD38.10gene.sig", t))
  ggsave(paste0(fn.fig, ".png"), w=2.5,h=7)
  ggsave(paste0(fn.fig, ".pdf"), w=2.5,h=7)
  
  # fn.fig = file.path(PROJECT_DIR, "figure_generation", sprintf("YF_Trial%d_CD38.10gene.sig_ROC",t))
  # ggsave(paste0(fn.fig,".png"), w=4,h=4)
  # ggsave(paste0(fn.fig,".pdf"), device=pdf, w=4,h=4)
  
  
  # fn.fig = file.path(PROJECT_DIR, "figure_generation", sprintf("YF_Trial%d_CD38.10gene.sig_BOXplot",t))
  # ggsave(paste0(fn.fig,".png"), w=4,h=4)
  # ggsave(paste0(fn.fig,".pdf"), device=pdf, w=4,h=4)
}
