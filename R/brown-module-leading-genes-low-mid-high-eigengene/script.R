library(WGCNA)
library(data.table)
library(pROC)

source(file.path(PROJECT_DIR, "R/brown-module-leading-genes-low-mid-high-eigengene/input.R"))

dn.out = file.path(PROJECT_DIR, "generated_data/brown-module-leading-genes-low-mid-high-eigengene")
dir.create(dn.out, showWarnings = F)


process.flu = function(tab, tab.info, fn, standardize) {

  m = t(as.matrix(tab[,-1,wi=F]))  # sample (row) by gene (column)
  colnames(m) <- toupper(tab[[1]])

  mask = rownames(m) %in% tab.info[which(Response %in% c("high","middle", "low"))]$SubjectID

  m = m[mask,]
  hi  = as.numeric((rownames(m) %in% tab.info[which(Response %in% "high")]$SubjectID))
  mid = as.numeric((rownames(m) %in% tab.info[which(Response %in% "middle")]$SubjectID))
  lo  = as.numeric((rownames(m) %in% tab.info[which(Response %in% "low")]$SubjectID))

  if(standardize) {
    mm = apply(m, 2, scale)
    rownames(mm) <- rownames(m)
    m  <- mm
  }
  
  mod.mask = (colnames(m) %in% brown.LE.isv050)

  mColorh = rep("grey", ncol(m))
  mColorh[mod.mask] <- "brown"
  
  MEList = moduleEigengenes(m, colors = mColorh)
  MEs = MEList$eigengenes
  MEs = as.data.table(MEs)
  MEs[,Sample:=rownames(m)]
  MEs[,Response:=0]
  MEs[which(mid== 1),Response:=1]
  MEs[which(hi == 1),Response:=2]

  fwrite(file=fn, MEs, quote=T)
}


fn = sprintf("%s/SDY212-zscore-AUC-low-mid-high-brown-leading-edge-wrt-CD38plusCellinFlu-%dgene.csv", dn.out, length(brown.LE.isv050))
process.flu(sdy212.d0, sdy212.samples, fn, standardize=T)

fn = sprintf("%s/SDY400-zscore-AUC-low-mid-high-brown-leading-edge-wrt-CD38plusCellinFlu-%dgene.csv", dn.out, length(brown.LE.isv050))
process.flu(sdy400.d0, sdy400.samples, fn, standardize=T)

fn = sprintf("%s/SDY404-zscore-AUC-low-mid-high-brown-leading-edge-wrt-CD38plusCellinFlu-%dgene.csv", dn.out, length(brown.LE.isv050))
process.flu(sdy404.d0, sdy404.samples, fn, standardize=T)

fn = sprintf("%s/H1N1-zscore-AUC-low-mid-high-brown-leading-edge-wrt-CD38plusCellinFlu-%dgene.csv", dn.out, length(brown.LE.isv050))
process.flu(H1N1.d0, H1N1.samples, fn, standardize=T)

