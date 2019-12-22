source(file.path(PROJECT_DIR, "R/functions/load_sig.r"))

fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_genes.txt")
fn.cd38.cor.all = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_all.genes.txt")
df.rob = read.table(fn.cd38.cor.all, sep="\t", header = TRUE, row.names = NULL)

fn.stab = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_genes_stability.txt")
df.stab = read.table(fn.stab, sep="\t", header = TRUE, row.names = 1)

# test signature
ng=20
gene.sig = load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=ng)

isv.seq = seq(0,1,by=0.05)
mat.rank = matrix(nrow=length(isv.seq), ncol=length(gene.sig))
rownames(mat.rank) = isv.seq
colnames(mat.rank) = gene.sig

df.comb = cbind(df.rob, df.stab) %>% 
  arrange(-cor.mean.sd.ratio)

for(i in seq_along(isv.seq)) {
  gi = df.comb$ISV >= isv.seq[i]
  mat.rank[i,] = match(gene.sig, df.comb$gene[gi]) / sum(gi)
}
isv.ng = sapply(isv.seq, function(x)sum(df.comb$ISV>=x))

df.rank = mat.rank %>% as.data.frame() %>% 
  tibble::rownames_to_column("ISV.th") %>% 
  gather("gene","rank",-ISV.th) %>% 
  mutate(ISV.th = as.numeric(ISV.th), gene = factor(gene, levels=gene.sig))

ggplot(df.rank, aes(ISV.th, rank, group=gene, col=gene)) +
  geom_vline(xintercept = 0.75, col="black", lty = 1) +
  geom_line(size=1) + 
  geom_point(fill="white", size=1, shape=21, stroke=1) +
  scale_x_reverse() +
  scale_y_reverse() +
  xlab("ISV threshold") +
  ylab("relative rank") +
  theme_bw() +
  theme(legend.key.size = unit(0.4, "cm"))

fn.fig = file.path(PROJECT_DIR, "figure_generation", "CD38.20gene.sig.rank.rel_vs_ISV.th")
ggsave(paste0(fn.fig,".png"), w=5,h=4)
ggsave(paste0(fn.fig,".pdf"), w=5,h=4)

