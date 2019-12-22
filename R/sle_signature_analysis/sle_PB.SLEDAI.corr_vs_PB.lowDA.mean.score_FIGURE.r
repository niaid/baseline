fn.pb = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_PB.DC_SLEDAI_corr_score_PG234_lowDA.txt")
fn.pb.low = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_lowDA_PB.DC_ge_sig_score_subjects.txt")

df.pb = fread(fn.pb)
df.pb.low = fread(fn.pb.low)

df = inner_join(df.pb, df.pb.low, by="SUBJECT") %>%
  dplyr::filter(PG %in% 2:4) %>% 
  mutate(PG = factor(PG)) %>% 
  mutate(PG.group = ifelse(PG %in% 2:4, PG.group.names[1], PG.group.names[2])) %>% 
  mutate(PG2 = ifelse(PG %in% 2:3, "2/3", ifelse(PG %in% 4, "4", "other"))) %>% 
  mutate(PG = factor(PG), PG2 = factor(PG2), PG.group = factor(PG.group, levels=PG.group.names))


cc = cor.test(df$PB_score_mean, df$PB_SLEDAI_corr_score, method="pearson")
df.cc = data.frame(label=sprintf("r = %.3f\np = %.2g",cc$estimate,cc$p.value),
                   x=-Inf, y=Inf)
df.cc2 = df %>% split(.$PG2) %>% 
  map(~cor.test(.$PB_score_mean, .$PB_SLEDAI_corr_score, method="pearson")) %>% 
  map_df(~data.frame(label=sprintf("r = %.3f\np = %.2g",.$estimate,.$p.value),
                     x=-Inf, y=Inf)) %>% 
  mutate(PG2 = levels(df$PG2) %>% factor())

PG2.clr = c("steelblue", "grey80")

p = ggplot(df, aes(PB_score_mean, PB_SLEDAI_corr_score)) +
  geom_point(shape=21, size=3, aes(fill=PG2), col="white") + #geom_smooth(aes(group=1),method="lm") +
  scale_fill_manual(values = PG2.clr, guide="none") +
  geom_text(data=df.cc, aes(label=label), x=-Inf, y=Inf, vjust=1.5, hjust=-0.1, col="black") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) + theme_bw() + theme(legend.key = element_blank()) + 
  xlab("Mean PB score at low DA") + 
  ylab("DA-associated\nchange in plasmablasts") +
  ggtitle("SLE patients with plasmablast\nsignature during flare") +
  theme_bw()
p
fn.fig = file.path(PROJECT_DIR, "figure_generation", "SLE_PB.DC.corr_vs_PB.mean.score_PG234_190628")
ggsave(paste0(fn.fig,".png"), w=2.6,h=3)
ggsave(paste0(fn.fig,'.pdf'), w=2.6,h=3)

