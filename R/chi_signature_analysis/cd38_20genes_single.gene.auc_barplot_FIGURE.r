library(pROC)

n_genes = 20

fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_genes.txt")
cc.rob = fread(fn.cd38.cor) %>% 
  arrange(desc(cor.mean.sd.ratio)) %>% 
  top_n(n_genes, cor.mean.sd.ratio) %>% 
  mutate(gene = factor(gene, levels=gene[order(.$cor.mean.sd.ratio)]))

# df.auc.one = data.frame(gene=factor(gene.sig, levels=rev(gene.sig)), AUC=auc.one)
ggplot(cc.rob, aes(x=gene, y=auc.gene)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept=0.5, lty=2, size=1, col="red") +
  geom_vline(xintercept=10.5, lty=2, size=1, col="black") +
  ylab("AUC") +
  coord_flip() + theme_bw()

fn.fig = file.path(PROJECT_DIR, "figure_generation", "CHI_CD38.20gene_single.gene.AUC")
ggsave(paste0(fn.fig,".png"), w=4, h=4)
ggsave(paste0(fn.fig,".pdf"), w=4, h=4)

ggplot(cc.rob, aes(x=gene, y=cor.mean.sd.ratio)) + 
  geom_bar(stat = "identity") +
  # geom_hline(yintercept=0.5, lty=2, size=1, col="red") +
  geom_vline(xintercept=10.5, lty=2, size=1, col="black") +
  ylab("Correlation coefficient, mean / SD") +
  coord_flip() + theme_bw()

fn.fig = file.path(PROJECT_DIR, "figure_generation", "CHI_CD38.20gene_single.gene.cor.mean.sd")
ggsave(paste0(fn.fig,".png"), w=4, h=4)
ggsave(paste0(fn.fig,".pdf"), w=4, h=4)


# plot all genes

cc.rob = fread(fn.cd38.cor) %>% 
  arrange(desc(cor.mean.sd.ratio)) %>% 
  # top_n(n_genes, cor.mean.sd.ratio) %>% 
  mutate(gene = factor(gene, levels=gene[order(.$cor.mean.sd.ratio)]))
cc.rob$sig = 0
cc.rob$sig[1:10] = 1

ggplot(cc.rob, aes(x=gene, y=cor.mean.sd.ratio, fill=factor(sig))) + 
  geom_col() +
  # geom_hline(yintercept=0.5, lty=2, size=1, col="red") +
  # geom_vline(xintercept=10.5, lty=2, size=1, col="black") +
  scale_fill_manual(values=c("black","red"), guide=F) +
  xlab("Genes ranked by correlation with CD20+CD38++ cell frequency") +
  ylab("") +
  # coord_flip() + 
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.grid = element_blank())

fn.fig = file.path(PROJECT_DIR, "figure_generation", "CHI_CD38.all.gene_single.gene.cor.mean.sd")
ggsave(paste0(fn.fig,".png"), w=5, h=1.5)
ggsave(paste0(fn.fig,".pdf"), w=5, h=1.5)
