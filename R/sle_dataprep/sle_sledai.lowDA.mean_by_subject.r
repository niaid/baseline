fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle_lowDA.txt")
info = fread(fn.si, data.table=F)

sledai = info %>% 
  dplyr::select(SUBJECT, SLEDAI) %>% 
  group_by(SUBJECT) %>% 
  summarise(SLEDAI.lowDA.mean = mean(SLEDAI), n.samples=n())

fn = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_SLEDAI.lowDA.mean.txt")
fwrite(sledai, fn, sep="\t")
