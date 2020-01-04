library(corrplot)

dn.fig = "figures/sig_ranked_correlation"
dir.create(dn.fig, showWarnings = F, recursive = T)

df.scores = fread("results/Score_table.txt")
names(df.scores) = sub("microarray", "PBMC", names(df.scores))

score.mat = df.scores %>% 
  dplyr::select(-response) %>% 
  tibble::column_to_rownames("subject") %>% data.matrix()

cc = cor(score.mat, method = "spearman", use = "pairwise.complete.obs")

cm.cor = pals::brewer.rdbu(10) %>% rev()

ccm = cor.mtest(score.mat, method = "spearman")
ccp = ccm$p
ord = c(13, 12, 1, 7, 14, 6, 4, 5, 2, 3, 9, 10, 8, 11)

fn.cor = glue::glue("{dn.fig}/Score_Spearman_corrplot_sig")
png(paste0(fn.cor, ".png"), w=500, h=600)
corrplot(cc[ord,ord], col=cm.cor,
         tl.col="black", cl.pos="b", order="original", 
         p.mat=ccp[ord,ord], insig="label_sig", sig.level=c(0.001, 0.01, 0.05), pch.col="white", pch.cex = 1)
dev.off()
pdf(paste0(fn.cor, ".pdf"), w=5, h=6, useDingbats = F)
corrplot(cc[ord,ord], col=cm.cor,
         tl.col="black", cl.pos="b", order="original", 
         p.mat=ccp[ord,ord], insig="label_sig", sig.level=c(0.001, 0.01, 0.05), pch.col="white", pch.cex = 1)
dev.off()

