fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info.txt")
info = fread(fn.si) %>% 
  mutate(subject = as.character(subject))

fn.demo = file.path(PROJECT_DIR, "data", "CHI", "phenotypes", "CHI_demographics.txt")
df.demo = fread(fn.demo) %>% 
  mutate(subject = as.character(subject))

fn.titer = file.path(PROJECT_DIR, "data", "CHI", "phenotypes", "titer_processed.txt")
df.titer = fread(fn.titer)  %>% 
  mutate(Subject = as.character(Subject)) %>% 
  mutate(Response = ifelse(adjMFC_class==0, "low",
                           ifelse(adjMFC_class==2, "high", "middle")))

info2 = info %>% 
  left_join(df.demo, by="subject") %>% 
  left_join(df.titer, by=c("subject"="Subject"))

fn.si.2 = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info_2.txt")
fwrite(info2, fn.si.2, sep="\t", quote=T)
