library(pROC)

# load gene expression data
fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix_gene.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()
fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info_2_CD38hi.txt")
info = fread(fn.si)

si = info$time==0 # select day 0 samples
dat = dat[,si]
info = info[si,]


niter = 500
set.seed(123)
cc.rob.rand = matrix(NA,nrow=nrow(dat),ncol=niter)
rownames(cc.rob.rand) = rownames(dat)
colnames(cc.rob.rand) = paste0("i",1:niter)

ns = 2
cmb = combn(1:ncol(dat), ns)

for (i in 1:niter) {
  cat(i,"")
  # calculate robust regression removing ns samples --------
  cc.rob = matrix(NA,nrow=nrow(dat),ncol=ncol(cmb))
  for (ic in 1:ncol(cmb)) {
    ridx = cmb[,ic]
    irand = sample(1:(ncol(dat)-2))
    cc.rob[,ic] = cor(t(dat[,-ridx]), info$CD38hi[-ridx][irand], method="spearman", use="pairwise.complete.obs")
  } # ic
  cc.rob.mean = apply(cc.rob,1,mean)
  cc.rob.sd = apply(cc.rob,1,sd)
  cc.rob.1cv = cc.rob.mean / cc.rob.sd

  cc.rob.rand[,i] = cc.rob.1cv
} # i

fn.rand = file.path(PROJECT_DIR, "generated_data", "CHI", 
                    sprintf("gene_sig_random_%d_from_all.txt",niter))
fwrite(cc.rob.rand, file=fn.rand, sep="\t")
