library(data.table)
library(latex2exp)
library(gplots)
library(GGally)
library(WGCNA)
allowWGCNAThreads()

dn = file.path(PROJECT_DIR, "R/WGCNA-modules-from-SLE-low-DA")
dn.out = file.path(PROJECT_DIR, "generated_data/WGCNA-modules-from-SLE-low-DA")
dir.create(dn.out, showWarnings = F)

source(file.path(dn, "input.R"))

datExpr0.dt = dat[T]
# dim(datExpr0.dt)
# datExpr0.dt[1:2,1:5]
datExpr0 = as.data.frame(t(datExpr0.dt[,-1,wi=F]))
names(datExpr0) <- datExpr0.dt$ProbeID
rownames(datExpr0) <- colnames(dat)[-1]
colnames(datExpr0) <- dat[[1]]
# dim(datExpr0)
# datExpr0[1:2,1:5]


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}



#=====================================================================================
#
#  WGCNA
#
#=====================================================================================


datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)

# datExpr[1:2,1:5]
fn.prefix= sprintf("%s/SLE-low-%dsbj-%dprobes", dn.out, nrow(datExpr), ncol(datExpr))
# fn.prefix


  # similar to above, except using "signed hybrid" instead of "signed"
  ############################################
  # pickSoftThreshold following the author's example 
  # https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/developingcortex/
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, corFnc="bicor", networkType="signed hybrid", verbose = 5, blockSize=6000)
  fn = sprintf("%s-pickSoftThreshold-signed-hybrid.rda", fn.prefix)
  # fn
  save(file=fn, sft)

  # Plot the results:
  if(1) {
    fn = sprintf("%s-soft-thresholding-powers-signed-hybrid.png", fn.prefix)
    png(fn, w=800, h=500)
      par(mfrow = c(1,2));
      cex1 = 0.9;
      # Scale-free topology fit index as a function of the soft-thresholding power
      plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed-hybrid R^2",type="n",
           ylim=c(-1,1),
           main = paste("Scale independence"));
      #abline(h=c(0,0.25,0.5,0.75), col='gray', lty='dashed')
      abline(h=seq(0,1,0.2), col='gray', lty='dashed')
      text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      labels=powers,cex=cex1,col="red");
      # this line corresponds to using an R^2 cut-off of h
      abline(h=0.80,col="red", lty='dotted')
      # Mean connectivity as a function of the soft-thresholding power
      plot(sft$fitIndices[,1], sft$fitIndices[,5],
           xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
           main = paste("Mean connectivity"), log="y")
      text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
      dev.off()
  }
  
  softPower = 4;
  fn = sprintf("%s-TOM-softPower%d-signed-hybrid.rda", fn.prefix, softPower)
  ## the following step needs at least 16GB and takes a long time to run!
  netData <- blockwiseModules(datExpr,maxBlockSize=20000,
                              networkType="signed hybrid",
                              power=softPower,
                              mergeCutHeight=0.15,
                              nThreads=10,
                              saveTOMFileBase = fn,
                              saveTOMs=TRUE,
                              corType="bicor",
                              minModuleSize=20,
                              pamStage=FALSE, ## The module size and other parameters do not matter right now - we re-cut the tree in the next step 
                              reassignThreshold = 1e-10,
                              verbose=3,
                              deepSplit=2) 
  fn = sprintf("%s-netData-saved.rda", fn.prefix)
  save(file=fn, netData)  # not sure if we need it later...
  
  TOMfile = sprintf("%s-TOM-softPower%d-signed-hybrid.rda-block.1.RData", fn.prefix, softPower)
  load(TOMfile)
  dissTOM = 1-TOM
  geneTree= hclust(dissTOM,method="average");

  # Plot the resulting clustering tree (dendrogram)
  fn = sprintf("%s-softPower%d-geneTree-signed-hybrid.pdf", fn.prefix, softPower)
  pdf(fn)
    plot(geneTree, xlab="", sub="", 
         main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
  dev.off()
  fn = sprintf("%s-softPower%d-geneTree-signed-hybrid.rda", fn.prefix, softPower)
  save(file=fn, geneTree)

  cut_and_merge <- function(minModSize, fn.prefix) {
    #minModSize: minimum number of genes in each module
    mColorh <- mLabelh <- colorLabels <- NULL
    dthresh <- 0.15 # MEs are no more than 0.85 correlated, if they are then the modules are merged and the ME is re-calculated
    ds <- 2 # deep split parameter to determine how finely to cut the tree
    tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
      minClusterSize = minModSize, cutHeight = 0.99999, ## Consider high maximum joining heights for dendrogram
      deepSplit = ds, distM = as.matrix(dissTOM))
    merged <- mergeCloseModules(exprData = datExpr,colors = tree$labels,
                                cutHeight = dthresh)
    mColorh <- labels2colors(merged$colors)
    mLabelh <- "DS=2,MMS=200,DCOR=0.15"
    fn = sprintf("%s-signed-hybrid-minModSize%d-TOMcuts.Rdata", fn.prefix, minModSize)
    save(geneTree,mColorh,mLabelh,file=fn)
    fn = sprintf("%s-geneTree-signed-hybrid-TOMcuts-minModSize%d.pdf", fn.prefix, minModSize)
    pdf(fn) 
        plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Defined Modules")
    dev.off()
    rev(sort(table(mColorh)))
    #system("printf press RET to continue; read a")
 
    ## see if we can merge some of the modules if their expression profiles are similar
    # Calculate eigengenes
    MEList = moduleEigengenes(datExpr, colors = mColorh)
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs);
    rownames(MEs) <- rownames(datExpr)
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average");
    # Plot the result
    #sizeGrWindow(7, 6)
    fn = sprintf("%s-METree-signed-hybrid-minModSize%d.pdf", fn.prefix, minModSize)
    pdf(fn)
        plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
    dev.off()
    ## merging threshold
    #MEDissThres = 0.25  ## corresponds to correlation 0.75
    ## Plot the cut line into the dendrogram
    #abline(h=dthresh, col = "red")
    fn = sprintf("%s-signed-hybrid-network-TOMcut-minModSize%d-eigen-genes.rda", fn.prefix, minModSize)
    save(file=fn, MEs, MEList, tree, merged, mColorh)
  }
  minModSize=20
  cut_and_merge(minModSize, fn.prefix)

  # export module membership
  fn = sprintf("%s-signed-hybrid-network-TOMcut-minModSize%d-eigen-genes.rda", fn.prefix, minModSize)
  attach(fn)
  ls(2)
  tt = rev(sort(table(get("mColorh",2))))
  #length(tt)
  fn = sprintf("%s-minModSize%d-%dmodules.pdf", fn.prefix, minModSize, length(tt)-1)
  pdf(fn)
       tt = rev(sort(table(get("mColorh",2))))[-1] # grey removed
       par(mar=c(5,8,5,5))
       main = sprintf("Number of Genes in the %d Modules", length(tt))
       barplot(tt, horiz=T, las=2, log="x", xlim=c(20,1200),
               xlab="#Genes in Module", main=main, col=names(tt))
       #text(tt[1]/1.1, 0.6, tt[1], pos=2, col="black") # for grey
       for(i in 1:length(tt)) {
         text(tt[i], i*1.2-0.6, tt[i], pos=4, col="black")
       }
  dev.off()

  gene.in.module = as.data.table(cbind(names(datExpr),get("mColorh",2)))
  setnames(gene.in.module,1:2,c("Symbol","Module"))
  #gene.in.module = merge(gene.in.module, annot, by="ProbeID")
  # dim(gene.in.module)
  setkey(gene.in.module, Module, Symbol)
  fn = sprintf("%s-gene.in.module.minModSize%d.signed-hybrid.txt", fn.prefix, minModSize)
  write.table(file=fn, gene.in.module, na="", quote=F, row.names=F, sep="\t")

  ## correlation between eigengene and genes in the module
  if(1) {
    attach(file.path(dn.out, "SLE-low-34sbj-9601probes-signed-hybrid-network-TOMcut-minModSize20-eigen-genes.rda"))
    # dim(get("MEs",2))
    # head(get("MEs",2))
    # dim(dat)  # gene x sample
    # head(dat,1)
    # head(gene.in.module)
    for(mod in unique(colnames(get("MEs",2)))) {
      g = gene.in.module[which(Module %in% gsub("^ME","",mod))]$Symbol
      d = dat[which(gene %in% g)]
      dd = as.matrix(d[,-1,wi=F])
      rownames(dd) = d[[1]]
      stopifnot(colnames(dd) == rownames(MEs))
      eigen.gene = MEs[,mod,drop=F]
      cor.res = as.data.table(cbind(t(cor(eigen.gene, t(dd))),
                               t(cor(eigen.gene, t(dd), method="spearman"))),keep.rownames=T)
      setnames(cor.res,"rn","Gene")
      setnames(cor.res,2,"cor.Pearson")
      setnames(cor.res,3,"cor.Spearman")
      fn = sprintf("%s/WGCNA-%s-gene-correlation.txt", dn.out, mod)
      write.table(file=fn, cor.res, sep="\t", quote=F, row.names=F)
    }
  }

# source("export.R")
# system("sh ./export.sh")
# setwd(cwd)
