library(data.table)
library(gplots)

dn = file.path(PROJECT_DIR, "R/WGCNA-eigengene-heatmap")
dn.out = file.path(PROJECT_DIR, "figure_generation/SLE-Sig")
dir.create(dn.out, showWarnings = F)
fn.prefix = file.path(dn.out, "SLE-low-34sbj-9601probes")

attach(file.path(PROJECT_DIR, "generated_data/WGCNA-modules-from-SLE-low-DA/SLE-low-34sbj-9601probes-signed-hybrid-network-TOMcut-minModSize20-eigen-genes.rda"))

n.mod = ncol(MEs)

# module X subject (labels are NOT colored)
## rows are modules

m = as.matrix(MEs)
colnames(m) <- gsub("ME","",colnames(m))
rownames(m) <- rownames(MEs)
m = t(m)
m = m[setdiff(rownames(m), "grey"),]

library(ComplexHeatmap)
library(circlize)
# cm = pals::brewer.brbg(20) #avoid pals package
cm = c("#543005","#714107","#8E530B","#A96C1E","#C28735","#D3AA5F","#E2C788","#EEDBAC","#F5EACD","#F5F1E7",
       "#E8F2F0","#D0ECE8","#B0E0D9","#8BD1C6","#64B9AE","#3C9C94","#1F827A","#036860","#005248","#003C30")
hm = Heatmap(m, name = "Module score", col = cm)
fn = sprintf("%s-eigengene-heatmap.pdf", fn.prefix)
pdf(fn, w=6, h=4)
draw(hm)
dev.off()

