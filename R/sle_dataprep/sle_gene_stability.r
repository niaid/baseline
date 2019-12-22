fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene_sle_lowDA.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle_lowDA.txt")
info = fread(fn.si, data.table=F)

# stability calculation
ss = matrix(NA,nrow(dat), 2, dimnames=list(rownames(dat),c("ISV","WSV")))
pv = rep(NA,nrow(dat))
names(pv) = rownames(dat)
for (i in 1:nrow(dat)) {
  tmp = info %>% 
  dplyr::select(SUBJECT) %>% 
  mutate(value=dat[i, ])
  fit = aov(value~SUBJECT, data=tmp)
  ss[i,] = summary(fit)[[1]]["Sum Sq"][[1]]
  pv[i] = summary(fit)[[1]]["Pr(>F)"][[1]][1]
}

pv.bh = p.adjust(pv, "BH")
pv.bonf = p.adjust(pv, "bonferroni")
ssn = ss / rowSums(ss)

fn.ssn = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_lowDA_genes_stability.txt")

ssn %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% 
  mutate(pv = pv, pv.BH=pv.bh, pv.Bonf=pv.bonf) %>% 
  fwrite(fn.ssn, sep="\t", col.names=T, row.names=F, quote=F)
