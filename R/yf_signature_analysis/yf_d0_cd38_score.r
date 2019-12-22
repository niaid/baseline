library(pROC)
source(file.path(PROJECT_DIR, "R/functions/get_score.r"))

# load CD38 signature genes
fn.sig = file.path(PROJECT_DIR, "generated_data", "signatures", "CD38_ge_sig.txt")
cd38.genes = fread(fn.sig, header = F) %>% unlist(use.names=F)

# load gene expression data
fn.ge = file.path(PROJECT_DIR, "generated_data", "YF", "YF_GE_matrix_gene_day0.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "YF", "YF_sample_info_day0.txt")
info = fread(fn.si)

gi = toupper(rownames(dat)) %in% toupper(cd38.genes)
sum(gi)

XX = rep(NA,nrow(info))

for (trial in c(1,2)) {
  lbl = paste0("T",trial)
  print(lbl)
  si = info$Trial==trial
  # ihl = info$Response[si] %in% c("low","high")

  X = get_score(dat[gi,si])

  XX[si] = X
}


df.out  = info %>% dplyr::select(subject, Trial, Response) %>% 
  mutate(CD38_score = XX)

fn.out = file.path(PROJECT_DIR, "generated_data", "YF", "YF_cd38_ge_sig_score.txt")
fwrite(df.out, fn.out, sep="\t", quote=T)
