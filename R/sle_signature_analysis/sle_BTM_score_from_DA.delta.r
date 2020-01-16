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

res.all = matrix(NA, nrow = nrow(list.mod$MODULES), ncol = ncol(dat2))
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
  mutate(Group = factor(Group, levels = unique(df.btm$Group))) %>% 
  group_by(Group) %>% 
  tally()

hm = Heatmap(mat.btm, name = "-log10(p.adj)", cluster_columns = F, cluster_rows = F, split=pg,
        col = colorRamp2(c(0, 3), c("black", "#F3EC18"))) + 
  Heatmap(pg, name="PG", col = cm, show_row_names = T, show_heatmap_legend = F)
# for(i in 1:7){
#   decorate_heatmap_body("-log10(p.adj)", slice = i, {
#     for(j in 1:(nrow(df.n)-1)){
#       grid.lines(rep(sum(df.n$n[1:j])/sum(df.n$n),2), c(0,1), gp = gpar(col = "white", lwd=3, lty = 1))
#     }
#     })
# }
fn.hm1 = file.path(dn.out, "SLE_BTM_PG")
png(paste0(fn.hm1, ".png"), w=800, h=800)
draw(hm)
dev.off()
pdf(paste0(fn.hm1, ".pdf"), w=9, h=10)
draw(hm)
dev.off()


# Compressed heatmap

DF2 = DF %>% 
  gather("ID", "p", -SUBJECT, -PG) %>% 
  left_join(df.btm, by="ID") %>% 
  mutate(Group = factor(Group, levels=rev(unique(Group)))) %>% 
  group_by(PG, Group) %>% 
  summarise(p.mean = mean(p, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(PG = paste0("PG",PG))

ggplot(DF2, aes(PG, Group, col=p.mean, size=p.mean)) +
  geom_point() +
  scale_color_distiller(palette = "YlOrRd", direction = 1, 
                        guide = guide_colorbar(title="mean -log(FDR)", title.position = "top", title.hjust = 0.5)) +
  scale_size_continuous(range = c(0, 25), guide = F) +
  scale_x_discrete(position = "top") +
  xlab("SLE Patient Groups") +
  ylab(NULL) +
  theme_bw() + theme(legend.position = "bottom")
fn.hm2 = file.path(dn.out, "SLE_BTM_PG_compressed2")
ggsave(paste0(fn.hm2, ".png"), w=7, h=5)
ggsave(paste0(fn.hm2, ".pdf"), w=7, h=5, useDingbats=F)
