library(data.table)

if(1) {
  fn = file.path(PROJECT_DIR, "generated_data/SLE/SLE_lowDA_genes_stability.txt")
  SLE.stability.based.on.log2 = fread(fn)
}

## dat
if(1) {
  fn = file.path(PROJECT_DIR, "generated_data/SLE/SLE_lowDA_PG234_ge.mean_matrix.txt")
  pv.BH.cutoff = 0.05
  dat = fread(fn)
  dat = dat[which(gene %in% SLE.stability.based.on.log2[which(pv.BH<=pv.BH.cutoff)][['gene']])]
  dim(dat)
  rm(pv.BH.cutoff)
}

## dat.scores
if(1) {
  fn = file.path(PROJECT_DIR, "generated_data/SLE/SLE_PB.DC_SLEDAI_corr_score.txt")
  dat.scores = fread(fn)
  # head(dat.scores)
  setnames(dat.scores, "SUBJECT", "Subject")
  dat.scores = dat.scores[which(Subject %in% colnames(dat))]
  dim(dat.scores)
  PB.score.cutoff = 1.2
  PB.score.name = "PB_SLEDAI_corr_score"
}

rm(fn)
