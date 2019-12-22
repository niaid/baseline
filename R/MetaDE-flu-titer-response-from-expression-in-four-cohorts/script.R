library(data.table)
library(MetaDE)
library(ggplot2)

#library(help="MetaDE")
#?MetaDE.ES

source(file.path(PROJECT_DIR, "R/MetaDE-flu-titer-response-from-expression-in-four-cohorts/input.R"))

x = list(list(exp1, label1$class), list(exp2, label2$class),
         list(exp3, label3$class), list(exp4, label4$class))

set.seed(20170706)
ind.res <- ind.cal.ES(x, paired=rep(FALSE,4), nperm=100)
names(ind.res)
res = MetaDE.ES(ind.res, meta.method="REM")

res.dt = as.data.table(lapply(res, function(x){x}))
res.dt[,gene:=names(res[[1]])]
setkey(res.dt, FDR.REM)
setkey(res.dt, zval)
res.dt = res.dt[nrow(res.dt):1]
fn = sprintf("%s/MetaDE-Effect-Size-from-four-datasets-%dgenes.csv", dn.out, nrow(res.dt))
fwrite(file=fn, res.dt, quote=T)


library(fgsea)

ranks = res.dt$zval
names(ranks) = res.dt$gene

mod.names = unique(WGCNA$Module)
mod.names = setdiff(mod.names, "grey")
mods = lapply(mod.names, function(x){WGCNA[which(Module %in% x)]$Symbol})
names(mods) <- mod.names

res.fgsea = fgsea(mods, ranks, nperm=100)

fn = sprintf("%s/enrichment-of-%dWGCNA-modules-in-4flu-datasets-%dcommon-genes.csv", dn.out, length(mod.names), nrow(res.dt))
fwrite(file=fn, res.fgsea, quote=T)

