source("R/functions/get_score.r")

fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle_lowDA.txt")
info = fread(fn.si, data.table=F)

fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene_sle_lowDA.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()


fn.sig = file.path(PROJECT_DIR, "generated_data", "signatures", "IFN26_ge_sig.txt")
gene.sig = fread(fn.sig, header = F) %>% unlist(use.names=F)

gi = toupper(rownames(dat)) %in% toupper(gene.sig)
sum(gi)

df.score = cbind(
                  dplyr::select(info, SUBJECT, VISIT, CUMULATIVE_TIME, SAMPLE_NAME), 
                  data.frame(score=get_score(dat[gi,]))
                )

fn.sig = file.path(PROJECT_DIR, "generated_data", "IFN26", "SLE_lowDA_IFN26_ge_sig_score.txt")
fwrite(df.score, fn.sig, sep="\t", quote=F)

df.score.subj = df.score %>%
  dplyr::select(SUBJECT, score) %>% 
  mutate(SUBJECT = factor(SUBJECT,
                          levels=unique(SUBJECT[order(as.numeric(sub("SLE-","",SUBJECT)))]))) %>% 
  group_by(SUBJECT) %>%
  dplyr::summarise(score_mean=mean(score, na.rm=T)) %>%
  ungroup()

fn.sig.subj = file.path(PROJECT_DIR, "generated_data", "IFN26", "SLE_lowDA_IFN26_ge_sig_score_subjects.txt")
fwrite(df.score.subj, fn.sig.subj, sep="\t", quote=F)
