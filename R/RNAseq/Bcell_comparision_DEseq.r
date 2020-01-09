# load libraries
# library("BiocParallel")
library("DESeq2")
library("RColorBrewer")
library("ggplot2")

# directories for results and figures
dn.out = file.path(PROJECT_DIR, "generated_data/RNAseq")
dir.create(dn.out, showWarnings = F, recursive = T)
dn.fig = file.path(PROJECT_DIR, "figure_generation/RNAseq")
dir.create(dn.fig, showWarnings = F, recursive = T)


## reading data from STAR counts output
counts = read.table(file = file.path(PROJECT_DIR, "data/RNAseq/counts.tab"),header = FALSE, row.names=1,sep =" ")
counts_run2 = read.table(file = file.path(PROJECT_DIR, "data/RNAseq/counts_run2.tab"),header = FALSE, row.names=1,sep =" ")
stopifnot(identical(rownames(counts), rownames(counts_run2)) & ncol(counts)==ncol(counts_run2))

countData = counts[5:22431,] + counts_run2[5:22431,]
countData = data.matrix(countData)

colData = read.csv(file = file.path(PROJECT_DIR, "data/RNAseq/colData.csv"), header = TRUE)
colnames(countData) = paste(colData$Donor, colData$Subset, sep=".")

# colSums(counts[5:22431,])
# colSums(counts_run2[5:22431,])
# colSums(countData)

# colSums(counts)
# colSums(countData)/10^6
# counts[1,]/colSums(counts) #percent unmapped
# counts[3,]/colSums(counts) #percent with no feature
# colSums(counts + counts_run2)/10^6

# removing a samples because of very low read counts due to low RNA/bad libraries 
countData = countData[,-4]
colData = colData[-4,]


#Create DEseq object and run DEseq2

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Donor + Subset)
dds = DESeq(dds)
vst = varianceStabilizingTransformation(dds)
# plotPCA(vst, intgroup = "Subset")
# plotPCA(vst, intgroup = "Donor")
# resultsNames(dds)

resCD19pCD20pvsCD19pCD20pCD38p=results(dds, contrast = c("Subset","CD19pCD20pCD38p","CD19pCD20p"))
summary(resCD19pCD20pvsCD19pCD20pCD38p)

write.csv(file = file.path(dn.out, "resCD19pCD20pvsCD19pCD20pCD38p.csv"), resCD19pCD20pvsCD19pCD20pCD38p, quote = FALSE)

vstData = assay(vst)
write.csv(file = file.path(dn.out, "VarianceStabilizedTransformCounts_CD38subsets.csv"), vstData)

CD38correlatedGenes = toupper(fread(file.path(PROJECT_DIR, "generated_data/signatures/CD38_ge_sig.txt"), header=F)[[1]])
CD38correlatedGenes = CD38correlatedGenes[CD38correlatedGenes %in% rownames(resCD19pCD20pvsCD19pCD20pCD38p)]

CD38correlatedGenes_resCD19pCD20pvsCD19pCD20pCD38p = as.data.frame(resCD19pCD20pvsCD19pCD20pCD38p[which(rownames(resCD19pCD20pvsCD19pCD20pCD38p) %in% CD38correlatedGenes),])
write.csv(file = file.path(dn.out, "CD38correlatedGenes_resCD19pCD20pvsCD19pCD20pCD38p.csv"), CD38correlatedGenes_resCD19pCD20pvsCD19pCD20pCD38p, quote = FALSE)

resCD19pCD20pvsCD19pCD20pCD38p_df = as.data.frame(resCD19pCD20pvsCD19pCD20pCD38p)
resCD19pCD20pvsCD19pCD20pCD38p_df_CD38correlatedGenes = resCD19pCD20pvsCD19pCD20pCD38p_df[CD38correlatedGenes,] %>% 
  tibble::rownames_to_column("gene")
resCD19pCD20pvsCD19pCD20pCD38p_df_sig = resCD19pCD20pvsCD19pCD20pCD38p_df[which(resCD19pCD20pvsCD19pCD20pCD38p_df$padj <= 0.01 & log(resCD19pCD20pvsCD19pCD20pCD38p$baseMean) > 1),]

sig_GOcellactivation = c("KLF6","MZB1","ADAM9","ACTB","XBP1","ACTG1","BCL2","IGLL5","LILRB1","CD84","SLAMF6","ADRB2","CD9","FCGR2A","CD38","CD44","TNFRSF13B","UBASH3B","LEF1","BANK1","CDK6")
resCD19pCD20pvsCD19pCD20pCD38p_df_sig_GOcellactivation = resCD19pCD20pvsCD19pCD20pCD38p_df[sig_GOcellactivation,] %>% 
  tibble::rownames_to_column("gene")

ggplot(data=resCD19pCD20pvsCD19pCD20pCD38p_df[which(log(resCD19pCD20pvsCD19pCD20pCD38p_df$baseMean) > 1),], aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(color = "red", alpha=0.4, size=1.75) +
  geom_point(data=resCD19pCD20pvsCD19pCD20pCD38p_df_CD38correlatedGenes, aes(x=log2FoldChange, y=-log10(padj)), size=2,colour="black") +
  geom_text(data=resCD19pCD20pvsCD19pCD20pCD38p_df_CD38correlatedGenes, aes(x=log2FoldChange, y=-log10(padj), label=gene, size=1.2),colour="black",hjust=0,vjust=1) +
  xlim(c(-7, 7)) +
  ylim(c(0, 9.5)) +
  xlab("log2 fold change") +
  ylab("-log10 FDR") +
  scale_colour_brewer(palette="Set1") +
  theme_bw() +
  theme(legend.position = "none")
fn.fig = file.path(dn.fig, "RNAseq_volcano")
ggsave(paste0(fn.fig, ".png"), w=6, h=5)


ggplot(data=resCD19pCD20pvsCD19pCD20pCD38p_df, aes(x=log(baseMean), y=log2FoldChange)) +
  geom_point(shape = 20, col="black", alpha=0.4, size=1.75) +
  ylim(c(-10,10)) +
  xlab("log(baseMean)") +
  ylab("log2 Fold Change") +
  geom_point(data=resCD19pCD20pvsCD19pCD20pCD38p_df_CD38correlatedGenes, aes(x=log(baseMean), y=log2FoldChange), size=2,colour="red") +
  geom_text(data=resCD19pCD20pvsCD19pCD20pCD38p_df_CD38correlatedGenes, aes(x=log(baseMean), y=log2FoldChange, label=gene, size=1.2),colour="red",hjust=0,vjust=1) +
  geom_point(data=resCD19pCD20pvsCD19pCD20pCD38p_df_sig_GOcellactivation, aes(x=log(baseMean), y=log2FoldChange), size=2,colour="deepskyblue") +
  geom_text(data=resCD19pCD20pvsCD19pCD20pCD38p_df_sig_GOcellactivation, aes(x=log(baseMean), y=log2FoldChange, label=gene, size=1.2),colour="deepskyblue",hjust=0,vjust=1) +
  theme_bw() +
  theme(legend.position = "none")
fn.fig = file.path(dn.fig, "RNAseq_MAplot")
ggsave(paste0(fn.fig, ".png"), w=6, h=6)
ggsave(paste0(fn.fig, ".pdf"), w=6, h=6)
