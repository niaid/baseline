library(pROC)
source(file.path(PROJECT_DIR, "R/functions/get_score.r"))
source(file.path(PROJECT_DIR, "R/functions/load_sig.r"))

# load CD38 signature genes
fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_genes.txt")
gene.sig = load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=10)

# load gene expression data
fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix_gene.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info_2_CD38hi.txt")
info = fread(fn.si)


roc.df = data.frame()

for (tp in c(0,-7,70)) {
  cat(tp,"")
  
  si = with(info, time == tp & !is.na(CD38hi))
  sum(si)
  tdat = dat[,si]
  tinfo = info[si,] %>% mutate(Response = factor(Response, levels=c("low","high")))
  
  day.pt = paste0("day ",tp)

  gi = toupper(rownames(tdat)) %in% toupper(gene.sig)
  sum(gi)
  
  ihl = tinfo$Response %in% c("low","high")
  
  # real signature
  X = get_score(tdat[gi,])[ihl]
  Y = tinfo$adjMFC_class[ihl]
  
  auc0 = roc(Y,X, direction="<", quiet = T)$auc
  
  # removing genes
  for (i in seq_along(gene.sig)) {
    isig = toupper(rownames(tdat)) %in% toupper(gene.sig[-i])
    X = get_score(tdat[isig,])[ihl]

    r = roc(Y,X, direction="<", quiet = T)
    r.df = data.frame(gene=factor(gene.sig[i], levels=gene.sig),
                      auc=as.numeric(r$auc)-auc0) %>% 
            mutate(day = case_when(
              tp == 0 ~ "Baseline 1 (day 0)",
              tp == -7 ~ "Baseline 2 (day -7)",
              tp == 70 ~ "Baseline 3 (day 70)"
            ))
    roc.df = rbind(roc.df, r.df)
  }
}

roc.df = roc.df %>% 
  # mutate(day = factor(day, levels=c("day 0","day -7","day 70")))
  mutate(day = factor(day))

auc.max = 0.2 #round(max(abs(range(roc.df.1g$auc))), digits = 1)

ggplot(roc.df, aes(gene, auc, fill=day)) + geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~day, nrow=1) +
  geom_hline(yintercept = 0, size=1, col="black") +
  scale_fill_manual(values=c("grey65", "grey80","grey50")) +
  xlim(rev(levels(roc.df$gene))) + 
  ylim(-auc.max, auc.max) + ylab("AUC change") +
  coord_flip() + theme_bw() +
  theme(legend.position="none", strip.background = element_blank(), 
        strip.text.x = element_text(size=12))

fn.fig = file.path(PROJECT_DIR, "figure_generation", "CD38.10gene.sig_AUC_remove_1gene")
ggsave(paste0(fn.fig,".png"), w=12, h=4.5)
ggsave(paste0(fn.fig,".pdf"), w=12, h=4.5)


