fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix_gene.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info_2.txt")
info = fread(fn.si)

si = info$time %in% c(-7,0,70)
dat = dat[,si]
info = info[si,] %>% mutate(subject = as.character(subject))

# stability calculation
ss = matrix(NA,nrow(dat), 2, dimnames=list(rownames(dat),c("ISV","WSV")))
pv = rep(NA,nrow(dat))
names(pv) = rownames(dat)
for (i in 1:nrow(dat)) {
  tmp = info %>% dplyr::select(subject) %>% mutate(value=dat[i, ])
  fit = aov(value~subject, data=tmp)
  ss[i,] = summary(fit)[[1]]["Sum Sq"][[1]]
  pv[i] = summary(fit)[[1]]["Pr(>F)"][[1]][1]
}

pv.bh = p.adjust(pv, "BH")
pv.bonf = p.adjust(pv, "bonferroni")
ssn = ss / rowSums(ss)

fn.ssn = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_genes_stability.txt")

ssn %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% 
  mutate(pv = pv, pv.BH=pv.bh, pv.Bonf=pv.bonf) %>% 
  fwrite(fn.ssn, sep="\t", col.names=T, row.names=F, quote=F)
