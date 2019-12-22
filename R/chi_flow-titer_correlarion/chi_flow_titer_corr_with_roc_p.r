library(pROC)
library(cowplot)

fn.poB = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_percent_of_B_filtered.txt")
flow.poB = fread(fn.poB) %>% tibble::column_to_rownames("sample") %>% 
  data.matrix()

fn.info = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_sample_info_filtered.txt")
flow.info = fread(fn.info) %>% 
  mutate(subject = as.character(subject))

fn.titer = file.path(PROJECT_DIR, "data", "CHI", "phenotypes", "titer_processed.txt")
df.titer = fread(fn.titer) %>% 
  mutate(Subject = as.character(Subject)) %>% 
  mutate(Response = ifelse(adjMFC_class==0, "low",
                           ifelse(adjMFC_class==2, "high", "middle")))

df = flow.info %>% add_column(CD38hi=flow.poB[,"Gate3"]) %>% 
  left_join(df.titer, by=c("subject"="Subject")) %>% 
  dplyr::filter(Response %in% c("low","high"),
         time %in% c(-7,0,70)) %>% 
  mutate(Response = factor(Response,levels=c("low","high"))) %>% 
  mutate(time.point = case_when(
                        time == 0 ~ "Baseline 1 (day 0)",
                        time == -7 ~ "Baseline 2 (day -7)",
                        time == 70 ~ "Baseline 3 (day 70)"
                      ) %>% factor())

# saveRDS(df, "flow_data_new.rds")


# generate Boxplots
df.w = df %>%
  dplyr::select(time.point, CD38hi, Response) %>%
  group_by(time.point) %>%
  dplyr::summarise(w.pv = ifelse(all(is.na(CD38hi)), NA,
    wilcox.test(CD38hi[Response=="high"], CD38hi[Response=="low"], exact=F, alternative = "greater")$p.value)) %>%
  mutate(lbl = sprintf("p = %.2g",w.pv)) %>%
  ungroup()

p1 = ggplot(df, aes(x=Response, y=CD38hi, group=Response)) +
  geom_boxplot(aes(fill=time.point), alpha=1, outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
  facet_wrap(~time.point, nrow=1) +
  scale_fill_manual(values=c("grey70", "grey90","grey50", "black", "white"), name="Response",
                    breaks=c(0,2), labels=c("low","high")) +
  geom_text(data = df.w, aes(label=lbl), x=1.5,y=Inf,vjust=1.1, hjust=0.5, size=4, inherit.aes = F) +
  xlab("Response") + ylab("CD38++ cells (% of alive B cells)") +
  coord_cartesian(xlim=c(0.8, 2.2)) +
  theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
        axis.ticks = element_blank())


# Generate ROC plots
roc.df = data.frame()
df.text = data.frame()
N_perm = 1000
for (tp in c(0,-7,70)) {
  r = roc(Response~CD38hi, data=dplyr::filter(df,time==tp), direction="<", quiet = T)
  set.seed(123)
  r.null = vector("numeric", N_perm)
  for(i in 1:N_perm) {
    # iperm = sample(nrow(df))
    r.null[i] = roc(Response~sample(CD38hi), data=dplyr::filter(df,time==tp), direction="<", quiet = T)$auc
  }
  r.p = sum(r.null > r$auc) / N_perm
  r.df = data.frame(day=paste0("day ",tp), Specificity = r$specificities, Sensitivity=r$sensitivities) %>%
    arrange(Sensitivity)
  roc.df = rbind(roc.df, r.df)
  
  df.text = rbind(df.text, 
                  data.frame(day=paste0("day ",tp), label=sprintf("AUC = %.2f\np = %.2g",r$auc,r.p)))
}
df.text$x = rep(0.25,3)
df.text$y = rep(0.1,3)
df.text = df.text %>% 
  mutate(day = factor(day, levels=levels(roc.df$day)))


p2 = list()
for(tm in levels(df.text$day)) {
  p2[[tm]] = ggplot(roc.df %>% dplyr::filter(day==tm), aes(x=Specificity, y=Sensitivity)) + 
    geom_line(size=1) +
    geom_abline(intercept = 1, slope = 1, lty=2) +
    scale_x_reverse() + 
    geom_text(data=df.text %>% dplyr::filter(day==tm), aes(x=x,y=y,label=label)) +
    coord_fixed() + 
    theme_bw() + theme(panel.grid = element_blank())
}

bottom_row <- do.call("plot_grid", c(p2, align = "h", nrow=1))
plot_grid(p1, bottom_row, ncol=1, rel_heights = c(2, 1))

fn.fig = file.path(PROJECT_DIR, "figure_generation", "CHI_flow_vs_respones_3time_p")
ggsave(paste0(fn.fig, ".png"), w=7,h=7)
ggsave(paste0(fn.fig, ".pdf"), w=7,h=7)
