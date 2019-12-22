fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info.txt")
info = read_tsv(fn.si)

fn.pg = file.path(PROJECT_DIR, "data", "SLE", "phenotypes", "SLE_SUBJECT_PG.txt")
df.pg = read_tsv(fn.pg)

info = info %>% 
  mutate(DA=ifelse(SLEDAI<3,"low",ifelse(SLEDAI>7,"high","mid"))) %>% 
  left_join(df.pg, by="SUBJECT")

fn.si.2 = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2.txt")
fwrite(info, fn.si.2, sep="\t", quote=T)
