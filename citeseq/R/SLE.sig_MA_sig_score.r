source("R/functions/get_score.r")

# load SLE.sig signature genes
sig.list = readRDS("sig/sig.list.RDS")
gene.sig = sig.list$SLE.sig

# load gene expression data
fn.ge = "../generated_data/CHI/CHI_GE_matrix_gene.txt"
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = "../generated_data/CHI/CHI_sample_info_2_CD38hi.txt"
info = fread(fn.si, data.table=F)

si = with(info, time == 0) # day 0 only
tdat = dat[,si]
tinfo = info[si,]

gi = toupper(rownames(tdat)) %in% toupper(gene.sig)
sum(gi)

ihl = tinfo$Response %in% c("low","high")
X = get_score(tdat[gi,])[ihl]

df = data.frame(subject=tinfo$subject[ihl], SLE.sig=X)

# output scores
fn.out = "results/sig_scores/scores_SLE.sig_MA.txt"
fwrite(df, fn.out, sep="\t")
