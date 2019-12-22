library(lme4)
source("R/functions/get_score.r")

fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle.txt")
info = fread(fn.si, data.table = F)

info = info %>% 
  mutate(DA=factor(DA, levels=c("low","mid","high")), 
         SUBJECT = factor(SUBJECT, levels=
                            unique(SUBJECT[order(as.numeric(sub("SLE-","",SUBJECT)))])))

treatment.fields = c("STEROID_IV_CATEGORY","CYCLOPHOSPHAMIDE_CATEGORY",
                     "ORAL_STEROIDS_CATEGORY","MYCOPHENOLATE_CATEGORY",
                     "HYDROXYCHLOROQUINE_CATEGORY")
for (f in treatment.fields) {
  info[,f] = factor(info[,f])
}


fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene_sle.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

pb.file = "PB_DC.M4.11_ge_sig.txt"
pb.label = "PB.DC"

fn.pg = file.path(PROJECT_DIR, "data", "SLE", "phenotypes", "SLE_SUBJECT_PG.txt")
df.pg = read_tsv(fn.pg)

subj.lowDA = info %>% 
               dplyr::select(SUBJECT, DA) %>% 
               dplyr::filter(DA=="low") %>% 
               distinct()

  fn.sig = file.path(PROJECT_DIR, "generated_data", "signatures", pb.file)
  pb.genes = fread(fn.sig, header = F) %>% unlist(use.names=F)
  
  gi = toupper(rownames(dat)) %in% toupper(pb.genes)
  sum(gi)
  info = mutate(info, PB=get_score(dat[gi,]))
  form = paste0("PB"," ~ 1 + ",
                paste(treatment.fields,sep="",collapse=" + "), " + (SLEDAI|SUBJECT)")
  m.lmer = lmer(as.formula(form), data=info)
  mre = ranef(m.lmer, condVar=T)
  x = mre$SUBJECT
  pv = attr(x, "postVar")
  se = unlist(lapply(1:ncol(x), function(i) sqrt(pv[i, i, ])))
  mre.df = data.frame(y=unlist(x),
                          se=se,
                          ci=1.96*se,
                          nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                          SUBJECT=factor(rep(rownames(x), ncol(x)), levels=rownames(x)),
                          group=gl(ncol(x), nrow(x), labels=names(x))) %>%
    dplyr::filter(group %in% c("SLEDAI")) %>%
    mutate(z = y/ci)
  mre.z = mre.df %>%
    dplyr::select(SUBJECT,group,z) %>%
    spread(group,z) %>% 
    dplyr::rename(PB_SLEDAI_corr_score=SLEDAI) 
  
  df.out = left_join(mre.z, df.pg, by="SUBJECT")
  fn.sig = file.path(PROJECT_DIR, "generated_data", "SLE", 
                     sprintf("SLE_%s_SLEDAI_corr_score.txt",pb.label))
  fwrite(df.out, fn.sig, sep="\t", quote=F)
  
  df.out = df.out %>% dplyr::filter(PG %in% 2:4, SUBJECT %in% subj.lowDA$SUBJECT)
  
  fn.sig.2 = file.path(PROJECT_DIR, "generated_data", "SLE", 
                     sprintf("SLE_%s_SLEDAI_corr_score_PG234_lowDA.txt",pb.label))
  fwrite(df.out, fn.sig.2, sep="\t", quote=F)
