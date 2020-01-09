library(fgsea)

CD38.genes = fread(file.path(PROJECT_DIR, "generated_data/CHI/robust_corr_all.genes.txt"))
gene.ISV = fread(file.path(PROJECT_DIR, "generated_data/CHI/CHI_genes_stability.txt"))
brown.L.E = fread(file.path(PROJECT_DIR, "generated_data/fgsea_with_wgcna_modules/brown-leading-edge-isv050-87genes.txt"))
##
## More information on results columns from DESEQ2
## http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
##
## [1] "mean of normalized counts for all samples"             
## [2] "log2 fold change (MLE): condition treated vs untreated"
## [3] "standard error: condition treated vs untreated"        
## [4] "Wald statistic: condition treated vs untreated"        
## [5] "Wald test p-value: condition treated vs untreated"     
## [6] "BH adjusted p-values"

dn.out = file.path(PROJECT_DIR, "generated_data/RNAseq")
if(!file.exists(dn.out)) stop("Please run the script R/RNAseq/Bcell_comparision_DEseq.r")
dn.fig = file.path(PROJECT_DIR, "figure_generation/RNAseq")
dir.create(dn.fig, showWarnings = F, recursive = T)

process <- function(fn, lists, tab1, main.text, y.shows.count, ylab) {
  # lists = brown.L.E
  tab1 = na.omit(tab1)
  setkey(tab1, stat)
  
  set.seed(2017110101)

  ### ranked by wald stat pvalue
  setkey(tab1, pvalue)
  ranks = tab1$pvalue
  names(ranks) = tab1$V1
  res.fgsea = fgsea(lists, ranks, nperm=100)
  fn.csv = paste0(file.path(dn.out, fn), ".csv")
  fwrite(file = fn.csv, res.fgsea, quote=T)

  fn.pdf = paste0(file.path(dn.fig, fn), ".pdf")
  pdf(fn.pdf, width=3, height=3.5)
         par(lwd=2)
         main = sprintf("%s (p=%0.2f)", main.text, res.fgsea$pval[1])
         xlab = "Gene rank based on\ndifferential expression"
         o = order(ranks)[which(names(ranks) %in% lists[[1]])]
         e = ecdf(o)
         if(y.shows.count) {
           y = e(o) * length(o)
         } else {
           y = e(o) 
         }
         plot(o, y, type='l', xlim=c(1,length(ranks)), 
              xlab=xlab, ylab=ylab, main=main, cex.main=0.7, cex.lab=0.7, las=1)
         if(y.shows.count) {
           lines(c(1,length(ranks)), c(1, length(o)), col="gray")
         } else {
           lines(c(0,length(ranks)), c(0, 1), col="gray")
         }
  dev.off()
}

tab1 = fread(file.path(dn.out, "resCD19pCD20pvsCD19pCD20pCD38p.csv"))

fn = "Brown-leading-edge-enrichment-in-CD19pCD20pvsCD19pCD20pCD38p-DE-genes"
process(fn, brown.L.E, tab1, main="Enrichment of\nSLE-Sig genes",
        y.shows.count=T, ylab = "Number of genes from SLE-Sig")

#########

process2 <- function(fn, lists, cols, labs, tab1, main.text, ylab) {
  ##lists = brown.L.E
  tab1 = na.omit(tab1)
  setkey(tab1, stat)
  
  set.seed(2017110101)

  ### ranked by wald stat pvalue
  setkey(tab1, pvalue)
  ranks = tab1$pvalue
  names(ranks) = tab1$V1
  res.fgsea = fgsea(lists, ranks, nperm=100)
  fn.csv = paste0(file.path(dn.out, fn), ".csv")
  fwrite(file = fn.csv, res.fgsea, quote=T)

  fn.pdf = paste0(file.path(dn.fig, fn), ".pdf")
  pdf(fn.pdf, width=5, height=4)
         par(lwd=2, xpd=T, mar=par()$mar+c(0,0,0,7))
         main = main.text
         xlab = "Gene rank based on\ndifferential expression"
         for(k1 in 1:length(lists)) {
           o = order(ranks)[which(names(ranks) %in% lists[[k1]])]
           e = ecdf(o)
           y = e(o) 
           if(k1 == 1) {
           plot(o, y, type='l', xlim=c(1,length(ranks)), col=cols[k1],
                xlab=xlab, ylab=ylab, main=main, cex.main=0.7, cex.lab=0.7, las=1, ylim=c(0,1))
           } else {
             lines(o, y, col=cols[k1])
           }
         }
         lines(c(0,length(ranks)), c(0, 1), col="gray")
         labs.text = sprintf("%s (p=%0.2f)", labs, res.fgsea$pval)
         cat(labs.text)
         legend(max(o)*1.1,1, labs.text, col=cols, lty=rep("solid", length(cols)), cex=0.8)
  dev.off()
}

if(1) {
  process3 = function(ISV.cutoff) {
    x = CD38.genes[which(gene %in% gene.ISV[which(ISV >= ISV.cutoff)]$gene)]
    setkey(x, cor.mean.sd.ratio)
    x1 = c(list(tail(x,10)$gene),
           list(tail(x,30)$gene),
           list(tail(x,50)$gene))
    cols = c("red", "gold", "blue")
    labs = c("top 10", "top 30", "top 50")

    tab1 = fread(file.path(dn.out, "resCD19pCD20pvsCD19pCD20pCD38p.csv"))
    fn = sprintf("ISV%d-CD38-top-enrichment-in-CD19pCD20pvsCD19pCD20pCD38p-DE-genes", ISV.cutoff*100)
    tab1 = tab1[which(V1 %in% x$gene)]
    process2(fn, x1, cols, labs, tab1, main=sprintf("Starting with genes\nwith TSM ≥ %.2f", ISV.cutoff),
            ylab = "Fraction of genes")
  }
  process3(ISV.cutoff=0.5)
  process3(ISV.cutoff=0.75)
}
