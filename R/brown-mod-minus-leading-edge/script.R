library(WGCNA)
library(data.table)
library(pROC)
library(limma)
library(ggplot2)
library(gplots)


source(file.path(PROJECT_DIR, "R/brown-mod-minus-leading-edge/input.R"))

dn.out = file.path(PROJECT_DIR, "generated_data/brown-mod-minus-leading-edge")
dir.create(dn.out, showWarnings = F)


# perform column-wise standardization
process.flu = function(tab, tab.info, fn) {

  m = t(as.matrix(tab[,-1,wi=F]))  # sample (row) by gene (column)
  colnames(m) <- toupper(tab[[1]])
  rn <- rownames(m)

  mod.mask = (colnames(m) %in% brown.minus.LE)

  m = apply(m,2,scale)
  rownames(m) <- rn
  
  hi.lo.mask = rownames(m) %in% tab.info[which(Response %in% c("high","low"))]$SubjectID

  m = m[hi.lo.mask,]
  hi  = as.numeric((rownames(m) %in% tab.info[which(Response %in% "high")]$SubjectID))
  lo  = as.numeric((rownames(m) %in% tab.info[which(Response %in% "low")]$SubjectID))

  mColorh = rep("grey", ncol(m))
  mColorh[mod.mask] <- "brown"
  
  MEList = moduleEigengenes(m, colors = mColorh)
  MEs = MEList$eigengenes
  MEs = as.data.table(MEs)
  MEs[,Sample:=rownames(m)]
  MEs[,HighResponse:=0]
  MEs[which(hi == 1),HighResponse:=1]
  
  fwrite(file=fn, MEs, quote=T)
  
  roc.res = roc(hi, MEs[[1]], direction="<", quiet = T)

}


fn = sprintf("%s/SDY212-AUC-brown-minus-leading-edge-wrt-CD38plusCellinFlu-%dgene.csv", dn.out, length(brown.minus.LE))
process.flu(sdy212.d0, sdy212.samples, fn)

fn = sprintf("%s/SDY400-AUC-brown-minus-leading-edge-wrt-CD38plusCellinFlu-%dgene.csv", dn.out, length(brown.minus.LE))
process.flu(sdy400.d0, sdy400.samples, fn)

fn = sprintf("%s/SDY404-AUC-brown-minus-leading-edge-wrt-CD38plusCellinFlu-%dgene.csv", dn.out, length(brown.minus.LE))
process.flu(sdy404.d0, sdy404.samples, fn)

fn = sprintf("%s/H1N1-AUC-brown-minus-leading-edge-wrt-CD38plusCellinFlu-%dgene.csv", dn.out, length(brown.minus.LE))
process.flu(H1N1.d0, H1N1.samples, fn)

