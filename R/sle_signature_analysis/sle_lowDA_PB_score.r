source("R/functions/get_score.r")

fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle_lowDA.txt")
info = fread(fn.si, data.table=F)

fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene_sle_lowDA.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

pb.file = "PB_DC.M4.11_ge_sig.txt"
pb.label = "PB.DC"

fn.sig = file.path(PROJECT_DIR, "generated_data", "signatures", pb.file)
pb.genes = fread(fn.sig, header = F) %>% unlist(use.names=F)

gi = toupper(rownames(dat)) %in% toupper(pb.genes)
sum(gi)

df.pb.score = cbind(
                dplyr::select(info, SUBJECT, VISIT, CUMULATIVE_TIME, SAMPLE_NAME),
                data.frame(PB_score=get_score(dat[gi,]))
              )

fn.sig = file.path(PROJECT_DIR, "generated_data", "SLE", 
                   sprintf("SLE_lowDA_%s_ge_sig_score.txt",pb.label))
fwrite(df.pb.score, fn.sig, sep="\t", quote=F)

df.pb.score.subj = df.pb.score %>%
  dplyr::select(SUBJECT, PB_score) %>% 
  mutate(SUBJECT = factor(SUBJECT,
                          levels=unique(SUBJECT[order(as.numeric(sub("SLE-","",SUBJECT)))]))) %>% 
  group_by(SUBJECT) %>%
  dplyr::summarise(PB_score_mean=mean(PB_score, na.rm=T)) %>%
  ungroup()

fn.sig.subj = file.path(PROJECT_DIR, "generated_data", "SLE",
                        sprintf("SLE_lowDA_%s_ge_sig_score_subjects.txt",pb.label))
fwrite(df.pb.score.subj, fn.sig.subj, sep="\t", quote=F)
