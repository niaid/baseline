fn.cd38 = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_lowDA_cd38_ge_sig_score_subjects.txt")
fn.pb = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_PB.DC_SLEDAI_corr_score.txt")

df.cd38 = fread(fn.cd38)
df.pb = fread(fn.pb)

PG.group.names = c("SLE patients with plasmablast\nsignature during flare", "Other SLE patients")

df = inner_join(df.cd38, df.pb, by="SUBJECT") %>% 
  dplyr::filter(!is.na(PG)) %>% 
  mutate(PG.group = ifelse(PG %in% 2:4, PG.group.names[1], PG.group.names[2])) %>% 
  mutate(PG2 = ifelse(PG %in% 2:3, "2/3", ifelse(PG %in% 4, "4", "other"))) %>% 
  mutate(PG = factor(PG), PG2 = factor(PG2), PG.group = factor(PG.group, levels=PG.group.names))

df.cc = df %>% split(.$PG.group) %>% 
  map(~cor.test(.$CD38_score_mean, .$PB_SLEDAI_corr_score, method="pearson")) %>% 
  map_df(~data.frame(label=sprintf("r = %.3f\np = %.2g",.$estimate,.$p.value),
                     x=-Inf, y=Inf)) %>% 
  mutate(PG.group = levels(df$PG.group) %>% factor())

df.cc2 = df %>% split(.$PG2) %>% 
  map(~cor.test(.$CD38_score_mean, .$PB_SLEDAI_corr_score, method="pearson")) %>% 
  map_df(~data.frame(label=sprintf("r = %.3f\np = %.2g",.$estimate,.$p.value),
                     x=-Inf, y=Inf)) %>% 
  mutate(PG2 = levels(df$PG2) %>% factor())

PG2.clr = c("steelblue", "grey80", "black")

p = ggplot(df, aes(CD38_score_mean, PB_SLEDAI_corr_score)) +
  geom_point(aes(fill=PG2), size=3, shape=21, col="white") + 
  scale_fill_manual(values = PG2.clr, guide="none") +
  facet_wrap(~PG.group, nrow=1) +
  geom_text(data=df.cc, aes(label=label), x=-Inf, y=Inf, vjust=1.5, hjust=-0.1, col="black") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) + theme_bw() + theme(legend.key = element_blank()) + 
  xlab("Mean Baseline score at low DA") + 
  ylab("DA-associated\nchange in plasmablasts") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(10,"mm") )
p

fn.fig = file.path(PROJECT_DIR, "figure_generation", "SLE_PB.DC.corr_vs_CD38.mean.score_2groups_190628p")
ggsave(paste0(fn.fig,".png"), w=6,h=3)
ggsave(paste0(fn.fig,'.pdf'), w=5,h=3)

