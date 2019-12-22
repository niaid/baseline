library(clinfun)
source(file.path(PROJECT_DIR, "R/functions/gg_color_hue.r"))

fn.si = file.path(PROJECT_DIR, "generated_data", "YF", "YF_cd38_ge_sig_score.txt")
info = fread(fn.si) %>% 
  dplyr::filter(Response %in% c("low","middle","high")) %>% 
  mutate(Response = factor(Response, levels=c("low","middle","high")),
         label = ifelse(Trial==1, "Yellow Fever", paste0("Yellow Fever, Trial ", Trial)))

df.text = data.frame()

for (trial in 1:2) {
  tinfo = info %>% dplyr::filter(Trial == trial)
  
  w.pv = jonckheere.test(tinfo$CD38_score, as.numeric(tinfo$Response), alternative="inc")$p.value
  
  df.text = rbind(df.text, 
                  data.frame(Trial=trial, #label.auc=sprintf("AUC = %.2f\np = %.2g", r$auc, r.p), 
                             label.w.pv=sprintf("p = %.2g",w.pv)))
}

df.text$x = 0.25
df.text$y = 0.1

df.text = df.text %>% 
  mutate(label = ifelse(Trial==1, "Yellow Fever", paste0("Yellow Fever, Trial ", Trial)))

clr = gg_color_hue(4) %>% rev()

for (t in 1:2) {
  df.text.2 = df.text %>% 
    dplyr::filter(Trial==t)
  
  ggplot(info %>% dplyr::filter(Trial==t), aes(x=Response, y=CD38_score, group=Response)) +
    geom_boxplot(aes(fill=label), alpha=0.5, outlier.colour = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
    facet_wrap(~label, nrow=1) +
    scale_fill_manual(values=c("black", "white", "grey35", clr[4]), name="Response",
                      breaks=c(0,1,2), labels=c("low","middle","high")) +
    geom_text(data = df.text.2, aes(label=label.w.pv), x=1.5, y=Inf,vjust=1.1, hjust=0.5, size=4, inherit.aes = F) +
    xlab("Response") + ylab("Baseline signature score") +
    coord_cartesian(xlim=c(0.8, 3.2)) +
    theme_bw() + theme(legend.position="none") +
    theme(panel.border = element_blank(), strip.background = element_blank(), 
          strip.text.x = element_text(size=12),
          panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
          axis.ticks = element_blank())
    
  fn.fig = file.path(PROJECT_DIR, "figure_generation", sprintf("YF_Trial%d_CD38.10gene.sig_with_mid_Jonckheere", t))
  ggsave(paste0(fn.fig, ".png"), w=2.5,h=4)
  ggsave(paste0(fn.fig, ".pdf"), w=2.5,h=4)
  
}
