library(lme4)
source("R/functions/get_score.r")

fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle.txt")
info = fread(fn.si, data.table = F)

info = info %>% 
  mutate(DA=factor(DA, levels=c("low","mid","high")), 
         DAn = as.numeric(DA),
         SUBJECT = factor(SUBJECT, levels=
                            unique(SUBJECT[order(as.numeric(sub("SLE-","",SUBJECT)))])))

info = info %>% 
  mutate(DA2 = ifelse(DAn>1, "high", "low")) %>% 
  group_by(SUBJECT) %>% 
  mutate(DAn.min = min(DAn, na.rm=t), DAn.max = max(DAn, na.rm=t)) %>% 
  ungroup()

si = info$DAn.min==1 & info$DAn.max>1 & (info$DAn==info$DAn.min | info$DAn==info$DAn.max)
sum(si)
subj = info$SUBJECT[si] %>% as.character() %>% unique()


fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene_sle.txt")
dat = fread(fn.ge, data.table = F) %>% 
  gather("SAMPLE_NAME", "value", -gene) %>% 
  inner_join(info[si,] %>% dplyr::select(SAMPLE_NAME, SUBJECT, PG, DA2), by="SAMPLE_NAME")

datm = dat %>% 
  group_by(gene, SUBJECT, PG, DA2) %>% 
  summarise(value = mean(value, na.rm=T)) %>% 
  ungroup()

dat.fc = datm %>% 
  spread("DA2","value") %>% 
  mutate(Delta = high-low)

dat2 = dat.fc %>% 
  dplyr::select(-PG, -high, -low) %>% 
  spread(SUBJECT, Delta) %>% 
  tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.btm = file.path(PROJECT_DIR, "data", "SLE", "sle_figure6A_BTM.txt")
df.btm = fread(fn.btm) %>% 
  mutate(ID = paste0("DC.",ID))

fn.subj = file.path(PROJECT_DIR, "data", "SLE", "sle_figure6A_subjects.txt")
df.subj = fread(fn.subj)

library(tmod)
data(tmod)
list.mod = tmod[df.btm$ID]


fn.pg = file.path(PROJECT_DIR, "data", "SLE", "phenotypes", "SLE_SUBJECT_PG.txt")
df.pg = read_tsv(fn.pg)

dn.out = file.path(PROJECT_DIR, "figure_generation")
dir.create(dn.out, showWarnings = F)

res.all = matrix(NA, nrow = length(list.mod), ncol = ncol(dat2))
rownames(res.all) = df.btm$ID
colnames(res.all) = colnames(dat2)

for(s in 1:ncol(dat2)) {
  res = tmodCERNOtest(rownames(dat2)[order(dat2[,s], decreasing = T)] , 
                       qval = 999, mset=list.mod, order.by = "none")
  res.all[,s] = res$adj.P.Val
}

isubj = match(df.subj$SUBJECT, colnames(res.all))
isubj = isubj[!is.na(isubj)]
res.all = res.all[,isubj]


df.all = -log10(res.all) %>% t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("SUBJECT")
  

DF = left_join(df.all, df.pg, by="SUBJECT")

library(ComplexHeatmap)
library(circlize)

mat.btm = DF %>% 
  dplyr::select(-PG) %>% 
  tibble::column_to_rownames("SUBJECT") %>% 
  data.matrix()
colnames(mat.btm) = with(df.btm, glue::glue("{sub('DC.','',ID)} - {Title}"))
pg = DF %>% 
  dplyr::select(SUBJECT, PG) %>% 
  mutate(PG = factor(PG)) %>% 
  tibble::column_to_rownames("SUBJECT") %>% 
  data.matrix()

cm = palette()[-8]
df.n = df.btm %>% 
  mutate(Group = fct_inorder(Group)) %>% 
  group_by(Group) %>% 
  tally()

Heatmap(mat.btm, name = "-log10(p.adj)", cluster_columns = F, cluster_rows = F, split=pg,
        col = colorRamp2(c(0, 3), c("black", "#F3EC18"))) + 
  Heatmap(pg, name="PG", col = cm, show_row_names = F, show_heatmap_legend = F)
for(i in 1:7){
  decorate_heatmap_body("-log10(p.adj)", slice = i, {
    for(j in 1:(nrow(df.n)-1)){
      grid.lines(rep(sum(df.n$n[1:j])/sum(df.n$n),2), c(0,1), gp = gpar(col = "white", lwd=3, lty = 1))
    }
    })
}
fn.hm1 = file.path(dn.out, "SLE_BTM_PG.png")
dev.copy(png, fn.hm1, w=800, h=800)
dev.off()
fn.hm1 = file.path(dn.out, "SLE_BTM_PG.pdf")
dev.copy(pdf, fn.hm1, w=9, h=10)
dev.off()


# Compressed heatmap
DF2 = DF %>% 
  gather("ID", "p", -SUBJECT, -PG) %>% 
  left_join(df.btm, by="ID") %>% 
  mutate(Group = factor(Group, levels=unique(Group))) %>% 
  group_by(PG, Group) %>% 
  summarise(p.mean = mean(p, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(PG = paste0("PG",PG)) %>% 
  spread(PG, p.mean)
DF2.mat = DF2 %>% 
  tibble::column_to_rownames("Group") %>%
  data.matrix()

library(corrplot)

DF2.mat[DF2.mat>15] = 15
# cm.cor = pals::brewer.rdylbu(12) %>% rev() #avoid pals package
cm.cor = c("#313695","#436FB1","#6BA2CB","#9CCDE2","#CCE9F2","#F0F9D8",
           "#FEF0A9","#FDCD7E","#FA9C58","#EE613D","#D22B26","#A50026")
fn.hm2 = file.path(dn.out, "SLE_BTM_PG_compressed")
png(paste0(fn.hm2, ".png"), w=300, h=200)
corrplot(DF2.mat, is.corr = F, col=cm.cor, cl.lim=c(0,15), cl.length = 4,
         tl.col="black", cl.pos="b", order="original")
dev.off()
pdf(paste0(fn.hm2, ".pdf"), w=3, h=2, useDingbats=F)
corrplot(DF2.mat, is.corr = F, col=cm.cor, cl.lim=c(0,15), cl.length = 4,
         tl.col="black", cl.pos="b", order="original")
dev.off()
