fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene_sle_lowDA.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle_lowDA.txt")
info = fread(fn.si)

info = info %>%
  mutate(SUBJECT = factor(SUBJECT, levels=
                            unique(SUBJECT[order(as.numeric(sub("SLE-","",SUBJECT)))])))

si = info$PG %in% c(5,6,7)
sum(si)
dat2 = dat[,si]
info2 = info[si,]

df = info2 %>% 
  dplyr::select(SUBJECT) %>% 
  bind_cols(as.data.frame(t(dat2))) %>% 
  gather("gene","value",-SUBJECT) %>% 
  group_by(SUBJECT, gene) %>% 
  dplyr::summarise(mean=mean(value)) %>% 
  ungroup()

df.mean = df %>% 
  spread("SUBJECT","mean")

fn.out = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_lowDA_PG567_ge.mean_matrix.txt")
fwrite(df.mean, fn.out, 
            sep="\t", row.names=F, col.names=T, quote=F)
