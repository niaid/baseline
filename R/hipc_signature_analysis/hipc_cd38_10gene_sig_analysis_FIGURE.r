library(pROC)
library(cowplot)
source(file.path(PROJECT_DIR, "R/functions/gg_color_hue.r"))

fn.si = file.path(PROJECT_DIR, "generated_data", "HIPC", "HIPC_cd38_ge_sig_score.txt")
info = fread(fn.si)%>% 
  dplyr::filter(Response %in% c("low","high")) %>% 
  mutate(Response = factor(Response, levels=c("low","high")))

sdy = data.frame(Study = c("SDY212", "SDY400", "SDY404"),
                 label = c("Stanford 2008", "Yale 2012", "Yale 2011"))

info = info %>% 
  inner_join(sdy, by="Study")

roc.df = data.frame()
df.text = data.frame()
N_perm = 1000

for (study in unique(info$Study)) {
  tinfo = info %>% dplyr::filter(Study == study)
  
  r = roc(Response ~ CD38_score, data=tinfo, direction="<", quiet = T)
  set.seed(123)
  r.null = vector("numeric", N_perm)
  for(i in 1:N_perm) {
    r.null[i] = roc(Response~sample(CD38_score), data=tinfo, direction="<", quiet = T)$auc
  }
  r.p = sum(r.null > r$auc) / N_perm

  r.df = data.frame(Specificity = r$specificities, Sensitivity=r$sensitivities, Study=study) %>%
    arrange(Sensitivity)
  
  roc.df = rbind(roc.df, r.df)
  X = tinfo$CD38_score
  Y = tinfo$Response
  w.pv = wilcox.test(X[Y=="high"], X[Y=="low"], exact=F, alternative="greater")$p.value
  
  df.text = rbind(df.text, 
                  data.frame(Study = study, label.auc=sprintf("AUC = %.2f\np = %.2g", r$auc, r.p),
                             label.w.pv=sprintf("p = %.2g", w.pv)))
}

roc.df = roc.df %>% 
  mutate(Study=factor(Study)) %>%
  arrange(Sensitivity, Specificity) %>% 
  inner_join(sdy, by="Study")

df.text$x = rep(0.25,3)
df.text$y = rep(0.1,3)
df.text = df.text %>% 
  inner_join(sdy, by="Study")

# Generate boxplots

clr = gg_color_hue(4) %>% rev()

p1 = ggplot(info, aes(x=Response, y=CD38_score, group=Response)) +
  geom_boxplot(aes(fill=label), alpha=0.5, outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
  facet_wrap(~label, nrow=1) +
  scale_fill_manual(values=c("black", "white", clr[1:3]), name="Response",
                    breaks=c(0,2), labels=c("low","high")) +
  geom_text(data = df.text, aes(label=label.w.pv), x=1.5, y=Inf,vjust=1.1, hjust=0.5, size=4, inherit.aes = F) +
  xlab("Response") + ylab("Baseline signature score") +
  coord_cartesian(xlim=c(0.8, 2.2)) +
  theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
        axis.ticks = element_blank())

p2 = list()
for(tm in levels(df.text$label)) {
  p2[[tm]] = ggplot(roc.df %>% dplyr::filter(label==tm), aes(x=Specificity, y=Sensitivity)) + 
    geom_line(size=1) +
    geom_abline(intercept = 1, slope = 1, lty=2) +
    scale_x_reverse() + 
    geom_text(data=df.text %>% dplyr::filter(label==tm), aes(x=x,y=y,label=label.auc)) +
    coord_fixed() + 
    theme_bw() + theme(panel.grid = element_blank())
}

bottom_row <- do.call("plot_grid", c(p2, align = "h", nrow=1))
plot_grid(p1, bottom_row, ncol=1, rel_heights = c(2, 1))

fn.fig = file.path(PROJECT_DIR, "figure_generation", "HIPC_CD38.10gene.sig_p")
ggsave(paste0(fn.fig, ".png"), w=7,h=7)
ggsave(paste0(fn.fig, ".pdf"), w=7,h=7)
