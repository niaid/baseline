source(file.path(PROJECT_DIR, "R/functions/get_score.r"))

# load CD38 signature genes
fn.sig = file.path(PROJECT_DIR, "generated_data", "signatures", "IFN26_ge_sig.txt")
gene.sig = fread(fn.sig, header = F) %>% unlist(use.names=F)

df.out = data.frame()

for (study in c("SDY212", "SDY400", "SDY404")) {

  # load gene expression data
  fn.ge = file.path(PROJECT_DIR, "generated_data", "HIPC",
                    paste0(study, "_GE_matrix_gene.txt"))
  dat = fread(fn.ge, data.table = F) %>% 
    tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
    data.matrix()
  
  fn.si = file.path(PROJECT_DIR, "generated_data", "HIPC",
                    paste0(study, "_sample_info.txt"))
  info = fread(fn.si)
  
  gi = toupper(rownames(dat)) %in% toupper(gene.sig)
  sum(gi)
  
  if(sum(gi)>1) {
    X = get_score(dat[gi,])
  } else {
    X = dat[gi,]
  }
  
  sdy.out  = info %>% 
    mutate(score = X) %>% 
    dplyr::filter(time=="d0") %>% 
    dplyr::select(subject, Study, Response, score) 
  
  df.out = rbind(df.out, sdy.out)
  
}

fn.out = file.path(PROJECT_DIR, "generated_data", "IFN26", "HIPC_IFN26_ge_sig_score.txt")
fwrite(df.out, fn.out, sep="\t", quote=T)
