fn.cd38 = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_lowDA_cd38_ge_sig_score_subjects.txt")
fn.pb = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_PB.DC_SLEDAI_corr_score.txt")

df.cd38 = fread(fn.cd38)
df.pb = fread(fn.pb)

PG.group.names = c("SLE patients with plasmablast\nsignature during flare", "Other SLE patients")

df = full_join(df.cd38, df.pb, by="SUBJECT") %>% 
  dplyr::filter(!is.na(PG)) %>%
  mutate(PG = ifelse(PG %in% 2:3, PG, "Other")) %>% 
  mutate(PG = factor(PG)) %>% 
  mutate(PG.group = ifelse(PG %in% 2:3, PG.group.names[1], PG.group.names[2])) %>% 
  mutate(PG.group = factor(PG.group, levels=PG.group.names))

cc = df %>% dplyr::filter(PG %in% 2:3) %>% 
  cor.test(~ CD38_score_mean + PB_SLEDAI_corr_score, data=., method="spearman")

p = ggplot(df %>% dplyr::filter(PG %in% 2:3), aes(CD38_score_mean, PB_SLEDAI_corr_score)) +
  geom_point(aes(fill=PB_SLEDAI_corr_score), size=3, shape=21, col="white") + 
  scale_fill_gradient2(low = "blue", mid="grey", high="red", name = "DaCP", midpoint = 0) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) + theme_bw() + theme(legend.key = element_blank()) + 
  xlab("Mean Baseline score at low DA") + 
  ylab("DaCP") +
  theme_bw()
p
fn.fig = file.path(PROJECT_DIR, "figure_generation", "SLE_DaCP_vs_TGSig_col.DaCP_PG23")
ggsave(paste0(fn.fig,".png"), w=4.5, h=3.2)
ggsave(paste0(fn.fig,".pdf"), w=4.5, h=3.2, useDingbats=F)


par.cur = par()

ths = df$PB_SLEDAI_corr_score[!is.na(df$PB_SLEDAI_corr_score) & !is.na(df$CD38_score_mean) & df$PG %in% 2:3] %>% unique() %>% sort()
ths = ths[1:(length(ths)-2)]
cc = rep(NA, length(ths))
pv = cc
nn = cc
meth = "pearson"
for(i in seq_along(ths)){
  ith = !is.na(df$PB_SLEDAI_corr_score) & df$PB_SLEDAI_corr_score > ths[i] & (df$PG %in% 2:3)

  nn[i] = sum(ith)
  if(sum(!is.na(df$PB_SLEDAI_corr_score[ith]))<3 | sum(!is.na(df$CD38_score_mean[ith]))<3) next()

  cc[i] = cor.test(df$PB_SLEDAI_corr_score[ith], df$CD38_score_mean[ith], method = meth, 
              use = "pairwise.complete.obs")$estimate
  pv[i] = cor.test(df$PB_SLEDAI_corr_score[ith], df$CD38_score_mean[ith], method = meth, 
              use = "pairwise.complete.obs")$p.value
}

df.cc = tibble(th = ths, `N subjects` = nn, `Pearson's r` = cc, `-log10(p)` = -log10(pv))
df.cc = df.cc %>% 
  gather("metric","value", -th) %>% 
  mutate(metric = factor(metric, levels=names(df.cc)[-1]))
df.vl = data.frame(metric = "-log10(p)", x = -log10(c(0.05, 0.1)))
ggplot(df.cc, aes(value, th)) +
  geom_point(size=3) +
  facet_wrap(~metric, nrow=1, scales = "free_x") +
  xlab("") + ylab("DaCP threshold (lower bound)") + 
  geom_vline(data=df.vl, aes(xintercept=x), col="red", lty=2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12),
        panel.spacing = unit(10,"mm") )
  
fn.fig = file.path(PROJECT_DIR, "figure_generation", "SLE_DaCP_threshold2_correlations_PG23.png")
ggsave(paste0(fn.fig,".png"), w=5, h=3.5)
ggsave(paste0(fn.fig,".pdf"), w=5, h=3.5, useDingbats=F)
