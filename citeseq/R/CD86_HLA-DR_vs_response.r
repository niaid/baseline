library(Seurat)
library(pROC)

dn.fig = "figures/adt_vs_response"
dir.create(dn.fig, showWarnings = F, recursive = T)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

df.subj = h1@meta.data %>% mutate(response = str_remove(adjmfc.time, "d0 ")) %>% 
  dplyr::select(subject=sampleid, response) %>% 
  mutate(subject = factor(subject)) %>% 
  distinct() %>% 
  arrange(subject)

clustering = 1
h1 = SetAllIdent(h1, id = glue::glue("K{clustering}"))

cl = "C9"
tobj = SubsetData(h1, ident.use = cl)
dat = tobj@assay$CITE@data
meta = tobj@meta.data %>% 
  dplyr::rename(subject = sampleid) %>% 
  mutate(subject = factor(subject, levels=df.subj$subject))

gi = rownames(dat) %in% c("CD86_PROT", "HLA-DR_PROT")
sum(gi)
dat2 = aggregate(t(as.matrix(dat[gi,,drop=F])), list(subject = meta$subject), mean, drop=F)
dat2[is.na(dat2)] = 0
mat = dat2 %>% 
  tibble::column_to_rownames("subject") %>% 
  data.matrix()

df = dat2 %>%left_join(df.subj, by="subject")

test.col.levels = c("low","high")
df.test = data.frame()
for(prot in colnames(mat)) {
  X = mat[,prot]
  Y = df$response %>% factor(levels=test.col.levels)
  r = roc(Y,X, direction="<", quiet = T)
  w.pv = wilcox.test(X[Y==test.col.levels[2]], X[Y==test.col.levels[1]], alternative = "greater", exact=F)$p.value
  mean.diff = mean(X[Y==test.col.levels[2]], na.rm=T) - mean(X[Y==test.col.levels[1]], na.rm=T)
  t.stat = ifelse(length(unique(X)) == 1, 0,
                  t.test(X[Y==test.col.levels[2]], X[Y==test.col.levels[1]], alternative = "two")$statistic)
  
  df.test = rbind(df.test, 
                  data.frame(protein = prot, auc=r$auc, w.pv, mean.diff, t.stat, N = length(tobj@cell.names)) %>% 
                    mutate(cluster = cl, clustering)
  )
  
}

fwrite(df.test, "results/CITEseq_CD80_CD86_HLA-DR_vs_Response.txt", sep="\t")

df.fig = df %>% 
  gather("protein","value",-c(subject, response)) %>% 
  mutate(protein = sub("_PROT", "", protein)) %>% 
  mutate(response = factor(response, levels = test.col.levels))
df.txt = df.test %>% 
  mutate(label = glue::glue("p = {format(w.pv, digits=2)}")) %>% 
  mutate(protein = sub("_PROT", "", protein))
ggplot(df.fig, aes(response, value, fill=response)) +
  geom_boxplot(fill = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  scale_fill_manual(values = c("white","black")) +
  geom_text(data=df.txt, aes(label=label), x = 1.5, y=Inf, hjust = 0.5, vjust = 1.5, 
            col="black", inherit.aes = F) +
  xlab("Response") +
  ylab("Background-corrected relative expression") +
  facet_wrap(~protein, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(file.path(dn.fig, "CITEseq_CD86_HLA-DR_vs_Response.png"), w=5, h=4)
ggsave(file.path(dn.fig, "CITEseq_CD86_HLA-DR_vs_Response.pdf"), w=5, h=4, useDingbats = F)
