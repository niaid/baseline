library(pROC)

# load gene expression data
fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix_gene.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info_2.txt")
info = fread(fn.si) %>% 
  mutate(subject = as.character(subject))

fn.poB = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_percent_of_B_filtered.txt")
flow.poB = fread(fn.poB) %>% tibble::column_to_rownames("sample") %>% 
  data.matrix()

fn.info = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_sample_info_filtered.txt")
flow.info = fread(fn.info) %>% 
  mutate(subject = as.character(subject))

flow.df = flow.info %>% dplyr::select(sample) %>% add_column(CD38hi = flow.poB[,"Gate3"])

info = info %>% left_join(flow.df, by="sample")

fn.si.new = sub(".txt", "_CD38hi.txt", fn.si)
fwrite(info, fn.si.new, sep="\t", quote=T)

si = with(info, time == 0 & Response %in% c("low","high") & !is.na(CD38hi))
sum(si)
dat = dat[,si]
info = info[si,] %>% mutate(Response = factor(Response, levels=c("low","high")))

# select temporary stable genes --------------------------
fn.stab = "generated_data/CHI/CHI_genes_stability.txt"
ge.stab = fread(fn.stab, sep="\t")
gidx = ge.stab$ISV >= 0.75
sum(gidx)
stopifnot(identical(ge.stab$gene, rownames(dat)))
dat.ts = dat[gidx,]


# calculate robust regression removing ns samples --------
ns = 2
cmb = combn(1:nrow(info), ns)
# cc.rob = matrix(NA,nrow=nrow(edat.gene),ncol=ncol(cmb))
cc.rob = matrix(NA,nrow=sum(gidx),ncol=ncol(cmb))
for (ic in 1:ncol(cmb)) {
  cat(ic," ")
  ridx = cmb[,ic]
  cc.rob[,ic] = cor(t(dat[gidx,-ridx]), info$CD38hi[-ridx], method="spearman", use="pairwise.complete.obs")
}
rownames(cc.rob) = rownames(dat)[gidx]
cc.rob.rank = apply(-cc.rob,2,rank)
ntop = 20
cc.rob.ntop = rowSums(cc.rob.rank<=ntop)

cc.rob.mean = apply(cc.rob,1,mean)
cc.rob.median = apply(cc.rob,1,median)
cc.rob.sd = apply(cc.rob,1,sd)
cc.rob.1cv = cc.rob.mean / cc.rob.sd
cc.rob.ntop.rank = rank(cc.rob.ntop)

cc.rob.1cv.ord = order(cc.rob.1cv,decreasing = T)
cc.rob.1cv.sort = sort(cc.rob.1cv,decreasing = T)
head(cc.rob.1cv.sort, 10)

# calculate AUC for each gene
auc.one=matrix(nrow=nrow(dat.ts), ncol=3)
for (k in 1:nrow(dat.ts)) {
  X = dat.ts[k,]
  Y = info$Response
  auc.one[k, c(2,1,3)] = ci.auc(Y,X,direction="<", quiet=T)
}
rownames(auc.one) = rownames(dat.ts)
colnames(auc.one) = c("auc.gene", "auc.ci95.min", "auc.ci95.max")

# output
df.out = data.frame(cor.ntop20 = cc.rob.ntop / ncol(cmb) *100, 
                    cor.mean = cc.rob.mean, 
                    cor.sd = cc.rob.sd, 
                    cor.median = cc.rob.median,
                    cor.mean.sd.ratio = cc.rob.1cv, 
                    as.data.frame(auc.one)) %>%
  tibble::rownames_to_column("gene")

fn.cor = file.path(PROJECT_DIR, "generated_data","CHI","robust_corr_genes.txt")
fwrite(df.out, fn.cor, sep="\t", quote=T)

ng = 30
df.top = df.out %>% 
  select(-cor.median) %>% 
  arrange(desc(cor.mean.sd.ratio)) %>% 
  top_n(ng, cor.mean.sd.ratio)

names(df.top) = c(
  "Gene",
  "% of iterations in top 20",
  "Mean of correlation coefficients",
  "SD of correlation coefficients",
  "Mean to SD ratio",
  "AUC for predicting titer response",
  "AUC, lower bound of 95% CI",
  "AUC, upper bound of 95% CI"
)


fn.cor = file.path(PROJECT_DIR, "generated_data","CHI",sprintf("robust_corr_genes_top%d.txt", ng))
fwrite(df.top, fn.cor, sep="\t", quote=T)
