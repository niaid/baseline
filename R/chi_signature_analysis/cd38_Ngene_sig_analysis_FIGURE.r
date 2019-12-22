library(pROC)
source(file.path(PROJECT_DIR, "R/functions/load_sig.r"))
source(file.path(PROJECT_DIR, "R/functions/get_score.r"))
fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_genes.txt")

# load gene expression data
fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix_gene.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info_2.txt")
info = fread(fn.si)

ng.max = 30
roc.df = data.frame()

for (tp in c(0,-7,70)) {
  cat(tp,"")
  
  day.pt = paste0("day ",tp)
  # test signature
  for(ng in 2:ng.max) {
    # load CD38 signature genes
    gene.sig = load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=ng)
    gi = toupper(rownames(dat)) %in% toupper(gene.sig)
    sum(gi)
    
    iday = info$time==tp
    ihl = info$adjMFC_class[iday] %in% c(0,2)
    X = get_score(dat[gi,iday])[ihl]
    Y = info$adjMFC_class[iday][ihl]
    
    auc.ng = ci.auc(Y,X, direction="<", quiet=T) %>% matrix(nrow=1)
    colnames(auc.ng) = c("auc.ci95.min", "AUC", "auc.ci95.max")
    r.df = data.frame(day=day.pt, ng, as.data.frame(auc.ng))
    
    roc.df = rbind(roc.df, r.df)
  }
}


roc.df = roc.df %>% 
  mutate(day = factor(day, levels=c("day 0","day -7","day 70"))) %>% 
  mutate(ngx = ifelse(day=="day -7", ng-0.2, ifelse(day=="day 70", ng+0.2, ng)))

ggplot(roc.df, aes(ng, AUC, group=day, col=day)) + 
  # geom_linerange(aes(x=ngx, ymin=auc.ci95.min, ymax=auc.ci95.max)) +
  # geom_ribbon(aes(ymin=auc.ci95.min, ymax=auc.ci95.max, fill=day), alpha=0.2, col=NA) +
  geom_line(size=1) +
  geom_point(size=1, shape=21, fill="white", stroke=1) +
  scale_color_manual(values=c("grey65", "grey80","grey50")) +
  xlab("No. of genes") + ylab("AUC") + theme_bw() + 
  theme(legend.key = element_blank())

fn.fig = file.path(PROJECT_DIR, "figure_generation", 
                   sprintf("CHI_CD38.2-%d.genes_AUC_ci95",ng.max))
ggsave(paste0(fn.fig,".png"), w=5, h=4)
ggsave(paste0(fn.fig,".pdf"), w=5, h=4)
