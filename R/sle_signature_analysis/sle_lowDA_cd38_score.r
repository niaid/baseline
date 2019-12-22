source("R/functions/get_score.r")

fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle_lowDA.txt")
info = fread(fn.si, data.table=F)

fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene_sle_lowDA.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()


fn.sig = file.path(PROJECT_DIR, "generated_data", "signatures", "CD38_ge_sig.txt")
cd38.genes = fread(fn.sig, header = F) %>% unlist(use.names=F)

gi = toupper(rownames(dat)) %in% toupper(cd38.genes)
sum(gi)

df.cd38.score = cbind(
                  dplyr::select(info, SUBJECT, VISIT, CUMULATIVE_TIME, SAMPLE_NAME), 
                  data.frame(CD38_score=get_score(dat[gi,]))
                )

fn.sig = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_lowDA_cd38_ge_sig_score.txt")
fwrite(df.cd38.score, fn.sig, sep="\t", quote=F)

df.cd38.score.subj = df.cd38.score %>%
  dplyr::select(SUBJECT, CD38_score) %>% 
  mutate(SUBJECT = factor(SUBJECT,
                          levels=unique(SUBJECT[order(as.numeric(sub("SLE-","",SUBJECT)))]))) %>% 
  group_by(SUBJECT) %>%
  dplyr::summarise(CD38_score_mean=mean(CD38_score, na.rm=T)) %>%
  ungroup()

fn.sig.subj = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_lowDA_cd38_ge_sig_score_subjects.txt")
fwrite(df.cd38.score.subj, fn.sig.subj, sep="\t", quote=F)
