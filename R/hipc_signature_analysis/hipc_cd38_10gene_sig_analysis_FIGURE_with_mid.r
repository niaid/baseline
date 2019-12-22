library(clinfun)
source(file.path(PROJECT_DIR, "R/functions/gg_color_hue.r"))

fn.si = file.path(PROJECT_DIR, "generated_data", "HIPC", "HIPC_cd38_ge_sig_score.txt")
info = fread(fn.si)%>% 
  dplyr::filter(Response %in% c("low","middle","high")) %>% 
  mutate(Response = factor(Response, levels=c("low","middle","high")))

sdy = data.frame(Study = c("SDY212", "SDY400", "SDY404"),
                 label = c("Stanford 2008", "Yale 2012", "Yale 2011"))

info = info %>% 
  inner_join(sdy, by="Study")

df.text = data.frame()

for (study in unique(info$Study)) {
  tinfo = info %>% dplyr::filter(Study == study)
  
  j.pv = jonckheere.test(tinfo$CD38_score, as.numeric(tinfo$Response), alternative="inc")$p.value
  
  df.text = rbind(df.text, 
                  data.frame(Study = study, #label.auc=sprintf("AUC = %.2f\np = %.2g", r$auc, r.p),
                             label.j.pv=sprintf("p = %.2g", j.pv)))
}

df.text$x = rep(0.25,3)
df.text$y = rep(0.1,3)
df.text = df.text %>% 
  inner_join(sdy, by="Study")


# Generate boxplots

clr = gg_color_hue(4) %>% rev()

ggplot(info, aes(x=Response, y=CD38_score, group=Response)) +
  geom_boxplot(aes(fill=label), alpha=0.5, outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
  facet_wrap(~label, nrow=1) +
  scale_fill_manual(values=c("black", "white", "grey35", clr[1:3]), name="Response",
                    breaks=c(0,1,2), labels=c("low","middle","high")) +
  geom_text(data = df.text, aes(label=label.j.pv), x=1.5, y=Inf,vjust=1.1, hjust=0.5, size=4, inherit.aes = F) +
  xlab("Response") + ylab("Baseline signature score") +
  coord_cartesian(xlim=c(0.8, 3.2)) +
  theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
        axis.ticks = element_blank())

fn.fig = file.path(PROJECT_DIR, "figure_generation", "HIPC_CD38.10gene.sig_with_mid_Jonckheere")
ggsave(paste0(fn.fig, ".png"), w=7,h=4)
ggsave(paste0(fn.fig, ".pdf"), w=7,h=4)
