require(data.table)
require(tmod) 
require(gplots)
require(RColorBrewer)
require(scales)
library(methods)
library(ggplot2)

dn.out = file.path(PROJECT_DIR, "generated_data/WGCNA-modules-from-SLE-low-DA-tmod")
dir.create(dn.out, showWarnings = F)
dn.fig = file.path(PROJECT_DIR, "figure_generation/SLE-Sig")
dir.create(dn.fig, showWarnings = F)

WGCNA = fread(file.path(PROJECT_DIR, "generated_data/WGCNA-modules-from-SLE-low-DA/SLE-low-34sbj-9601probes-gene.in.module.minModSize20.signed-hybrid.txt"))
setkey(WGCNA, Module, Symbol)

mods = unique(WGCNA$Module)
res = as.list(rep(NA,length(mods)))
for(i in 1:length(mods)) {
  mod = mods[i]
  r = as.data.table(tmodHGtest(fg=WGCNA[mod]$Symbol, bg=WGCNA$Symbol))
  if(nrow(r)>0) {
    r[,WGCNA:=mod]
  }
  res[[i]] <- r
}

idx = which(sapply(res, nrow)>0)
res = rbindlist(res[idx])

res[,genes.in.overlap:=""]
res[,genes.in.bg:=""]
c1 = which(colnames(res) %in% "genes.in.overlap")
c2 = which(colnames(res) %in% "genes.in.bg")
for(i in 1:nrow(res)) {
  mod = res[i]$WGCNA
  g1 = getGenes(res$ID[i], fg=WGCNA[mod]$Symbol)$fg
  g2 = getGenes(res$ID[i], fg=WGCNA$Symbol)$fg
  set(res,i=i,j=c1,g1)
  set(res,i=i,j=c2,g2)
}

setnames(res,1,"BTM")  # avoid using ID as the first header to make excel happy

fn = file.path(dn.out, "WGCNA-BTM-enrichment-p.adjust.txt")
write.table(file=fn, res, row.names=F, quote=F, sep="\t")

res = res[!grepl("TBA", Title)]
res[,Desc:= sprintf("%s (%s)", BTM, Title)]
ID1 = unique(res$Desc)
ID2 = unique(res$WGCNA)
N1 = length(ID1)
N2 = length(ID2)

# ----- Barplot of brow module enrichment ----------
res.brown = res[WGCNA=="brown"] %>% 
  as.data.frame() %>% 
  arrange(adj.P.Val) %>% 
  mutate(Desc = factor(Desc, levels = rev(unique(Desc))))
ggplot(res.brown, aes(Desc, -log10(adj.P.Val))) +
  geom_col() +
  geom_hline(yintercept = -log10(0.01), col="red", lty=2) +
  xlab("") +
  ylab("-log10(FDR)") +
  coord_flip() +
  theme_bw()
ggsave(file.path(dn.fig, "BTM_enrichment_in_brown_module.pdf"), w=7, h=3)  

#--------------------------------------------------


m = matrix(rep(-1,N1*N2),nrow=N1, ncol=N2)
rownames(m) <- ID1
colnames(m) <- ID2

cn = "adj.P.Val"
for(k1 in 1:nrow(res)) {
  n1 = res[k1]$Desc
  n2 = res[k1]$WGCNA
  xx = -log10(res[k1][[cn]])
  m[n1,n2] = xx
}
m[m<0] = 0

rc = rep("black",nrow(m))
rc[m[,"brown"]>0] = "brown"

fn = file.path(dn.fig, "WGCNA-BTM-enrichment-p.adjust.pdf")
m.max = 10
pdf(fn, width=10, height=15)
    breaks = seq(min(m), m.max, length=51)
    hm = heatmap.2(m, Colv = F, trace='none', margin=c(10,20), breaks=breaks, 
              col=c("white",colorRampPalette(c("white","darkred"))(49)),
              ColSideColors=colnames(m),
              colRow=rc,
              dendrogram = "none",
              keysize = 0.75,
              lhei=c(1,6), lwid=c(2,6),
              key.xlab = "-log2(FDR)", key.title = "",
              density.info = "none")
dev.off()
