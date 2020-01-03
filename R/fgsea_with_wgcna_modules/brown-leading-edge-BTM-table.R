library(data.table)
library(tmod)
library(fgsea)
library(methods)

dn.out = file.path(PROJECT_DIR, "generated_data/fgsea_with_wgcna_modules/")
dir.create(dn.out, showWarnings = F)


fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_all.genes.txt")
df.cd38.cor = fread(fn.cd38.cor)
ranked = df.cd38.cor[, .(Gene=gene, cor.mean.sd.ratio)]
ranked[,Gene:=toupper(Gene)]

mods = fread(file.path(PROJECT_DIR, "generated_data/WGCNA-modules-from-SLE-low-DA/SLE-low-34sbj-9601probes-gene.in.module.minModSize20.signed-hybrid.txt"))
mod = toupper(mods[which(Module %in% "brown")]$Symbol)

r1 = 1:nrow(ranked)
r2 = r1[which(ranked$Gene %in% mod)]

m75   = strsplit(getGenes("LI.M75")$Genes, ",")[[1]]
m150  = strsplit(getGenes("LI.M150")$Genes, ",")[[1]]
m165  = strsplit(getGenes("LI.M165")$Genes, ",")[[1]]

if(1) {
  convert.to.dt = function(lst, name) {
    tmp = as.data.table(list(lst, rep(1, length(lst))))
    setnames(tmp,1, "Gene")
    setnames(tmp, "V2", name)
    return(tmp)
  }
  tab = Reduce(function(x,y){merge(x,y,all=T, by="Gene")},
               list(convert.to.dt(m75,"LI.M75"),   convert.to.dt(m150, "LI.M150"),
                    convert.to.dt(m165, "LI.M165")), convert.to.dt(mod, "brown"))
  fn.out = file.path(dn.out, "brown-mod-m75-m150-m165-genes.csv")
  fwrite(file=fn.out, tab, quote=T)
  tab = tab[which(!is.na(brown) & (!is.na(LI.M75) | !is.na(LI.M150) | !is.na(LI.M165)))]
  fwrite(file=sub(".csv", "-short.csv", fn.out), tab, quote=T)
}
