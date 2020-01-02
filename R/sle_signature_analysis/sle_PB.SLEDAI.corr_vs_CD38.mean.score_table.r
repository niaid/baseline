fn.cd38 = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_lowDA_cd38_ge_sig_score_subjects.txt")
fn.pb = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_PB.DC_SLEDAI_corr_score.txt")

df.cd38 = fread(fn.cd38)
df.pb = fread(fn.pb)

fn.pg = file.path(PROJECT_DIR, "data", "SLE", "phenotypes", "SLE_SUBJECT_PG.txt")
df.pg = read_tsv(fn.pg)

PG.group.names = c("SLE patients with plasmablast signature during flare", "Other SLE patients")

df = inner_join(df.cd38, df.pb, by="SUBJECT") %>% 
  dplyr::filter(!is.na(PG)) %>% 
  mutate(PG.group = ifelse(PG %in% 2:4, PG.group.names[1], PG.group.names[2])) %>% 
  mutate(PG = factor(PG), PG.group = factor(PG.group, levels=PG.group.names))

df.cc.s = df %>% split(.$PG.group) %>% 
  map(~cor.test(.$CD38_score_mean, .$PB_SLEDAI_corr_score, method="spearman")) %>% 
  map_df(~data.frame(Spearman.rho=.$estimate, Spearman.p = .$p.value)) %>% 
  mutate(PG.group = levels(df$PG.group) %>% factor())

df.cc.p = df %>% split(.$PG.group) %>% 
  map(~cor.test(.$CD38_score_mean, .$PB_SLEDAI_corr_score, method="pearson")) %>% 
  map_df(~data.frame(Pearson.cor=.$estimate, Pearson.p = .$p.value)) %>% 
  mutate(PG.group = levels(df$PG.group) %>% factor())

df.cc = inner_join(df.cc.p[,c(3,1,2)], df.cc.s, by="PG.group")
fn.out = file.path(PROJECT_DIR, "generated_data/SLE/SLE_PB.DC.corr_vs_CD38.mean.score_2groups.txt")
fwrite(df.cc, fn.out, sep="\t")
