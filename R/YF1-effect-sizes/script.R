library(data.table)
library(ggplot2)
library(effsize)
library(fgsea)

source(file.path(PROJECT_DIR, "R/YF1-effect-sizes/input.R"))

process <- function(fn, expr, annot) {
  response.grp0 = (colnames(expr) %in% annot[which(class==0)]$geo_accession)
  response.grp1 = (colnames(expr) %in% annot[which(class==1)]$geo_accession)
  
  es = apply(expr, 1, function(x) {
               x0 = x[response.grp0]
               x1 = x[response.grp1]
               r = cohen.d(x1, x0, hedges.correction=T)
               c(r$estimate, r$var, r$conf.int)
  })
  
  es = t(es)
  es = as.data.table(es, keep.rownames=T)
  setnames(es, 1, "Gene")
  setnames(es, 2, "Hedges.g")
  setnames(es, 3, "Var")
  setnames(es, 4, "Conf95.lower")
  setnames(es, 5, "Conf95.upper")

  es[,Gene:=toupper(Gene)]
  setkey(es, Gene)
  
  CD38.sig = toupper(fread(file.path(PROJECT_DIR, "generated_data/signatures/CD38_ge_sig.txt"), header=F)[[1]])

  fwrite(file=fn, es, quote=T)
  
  fn = paste0(gsub(".csv$", "", fn), "-CD38sigOnly.csv")
  fwrite(file=fn, es[CD38.sig,], quote=T)
}

fn = file.path(dn.out, "YF1-hedges-g.csv")
process(fn, exp1, annot1)

# fn = file.path(dn.out, "YF2-hedges-g.csv")
# process(fn, exp2, annot2)

## 
YF1 = fread(file.path(dn.out, "YF1-hedges-g.csv"))
YF1 = na.omit(YF1)
setkey(YF1, Hedges.g)
YF1 = YF1[rev(1:nrow(YF1))]
ranks = YF1$Hedges.g
names(ranks) = YF1$Gene

mods = fread(file.path(PROJECT_DIR, "generated_data/WGCNA-modules-from-SLE-low-DA/SLE-low-34sbj-9601probes-gene.in.module.minModSize20.signed-hybrid.txt"))
brown.mod = mods[which(Module %in% "brown")][,1,wi=F]
brown.mod = brown.mod[which(Symbol %in% names(ranks))]
# table(brown.mod[[1]] %in% names(ranks))
res.fgsea = fgsea(brown.mod, ranks, nperm=100)

fwrite(file=file.path(dn.out, "brown-mod-vs-YF1-genes-fgsea.csv"), res.fgsea, quote=T)

