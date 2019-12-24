library(lme4)

dn = file.path(PROJECT_DIR, "R/PB-score-vs-WGCNA-modules")
dn.out = file.path(PROJECT_DIR, "generated_data/PB-score-vs-WGCNA-modules")
dir.create(dn.out, showWarnings = F)
dn.fig = file.path(PROJECT_DIR, "figure_generation/SLE-Sig")
dir.create(dn.fig, showWarnings = F)

# load WGCNA modules data

fn = file.path(PROJECT_DIR, "generated_data/WGCNA-modules-from-SLE-low-DA/SLE-low-34sbj-9601probes-signed-hybrid-network-TOMcut-minModSize20-eigen-genes.rda")
load(fn)
MEs.tab = as.data.table(MEs, keep.rownames=T)
setnames(MEs.tab, "rn", "Subject") 

fn = file.path(PROJECT_DIR, "generated_data/SLE/SLE_PB.DC_SLEDAI_corr_score_PG234_lowDA.txt")
PB.score.description = "PB score"
PB.score.name = "PB_SLEDAI_corr_score"
PB.score.cutoff = 1.1 

# load scores

tab = fread(fn)
setnames(tab, "SUBJECT", "Subject")
tab = tab[which(Subject %in% MEs.tab$Subject)]

# merge WGCNA and TGSig scores
tab = merge(tab, MEs.tab, by="Subject") 

top.n = table(tab[[PB.score.name]] > PB.score.cutoff)[2]
tab[,PG:=NULL]

tab.PB1 = tab
setkey(tab.PB1, Subject)

end_point_name = PB.score.name 

dat = tab.PB1[T]
dat[,MEgrey:=NULL] 
setnames(dat, PB.score.name, "y")

dat[,Subject:=NULL]

cc1 = cor(dat$MEbrown, dat$y)
cc2 = cor(dat$MEbrown, dat$y, method="spearman")

# Correlation plot for brown module ---

main = sprintf("brown module vs %s\n(Pearson corr=%0.3f, Spearman corr=%0.3f)", PB.score.description,
               cc1, cc2)
xlab = PB.score.description 
ylab = "Brown module score"

pdf(file.path(dn.fig, "PB-score_vs_brown-module.pdf"), width=3, height=3.5)
plot(MEbrown ~ y, xlab=xlab, ylab=ylab, 
     cex.main = 0.7,
     main=main, pch=19, col="brown", data=dat, las=2)
fit = coef(lm(MEbrown ~ y, data=dat))
abline(h = seq(-0.2,0.4,0.2), col="gray", lty="dotted")
abline(v = seq(-1,3,1), col="gray", lty="dotted")
abline(fit[1], fit[2], col="black") 
dev.off()


# Empirical correlation ---

if(1) {
  N = 100000
  fn.cache = file.path(dn.out, "ecdf-shuffled-spearman-and-pearson-corrs.rda")
  if(file.exists(fn.cache)) {
    load(fn.cache)
  } else {
    pearson.shuffled = lapply(colnames(dat)[-1], function(mod) {
      set.seed(0)
      replicate(N,cor(sample(dat$y), dat[[mod]], method="pearson"))
    })
    pearson.ecdfs = lapply(pearson.shuffled, function(x) { ecdf(x) })
    spearman.shuffled = lapply(colnames(dat)[-1], function(mod) {
      set.seed(0)
      replicate(N,cor(sample(dat$y), dat[[mod]], method="spearman"))
    })
    spearman.ecdfs = lapply(spearman.shuffled, function(x) { ecdf(x) })
    names(pearson.shuffled) = colnames(dat)[-1]
    names(spearman.shuffled) = colnames(dat)[-1]
    names(pearson.ecdfs) = colnames(dat)[-1]
    names(spearman.ecdfs) = colnames(dat)[-1]
    save(file=fn.cache, pearson.ecdfs, spearman.ecdfs, pearson.shuffled, spearman.shuffled)
  }
  

  cc.pearson = cor(dat$y, dat[,-1], method="pearson")[1,]
  cc.spearman = cor(dat$y, dat[,-1], method="spearman")[1,]
  pvals.pearson = rep(NA, ncol(dat)-1)
  names(pvals.pearson) = colnames(dat)[-1]
  pvals.spearman = rep(NA, ncol(dat)-1)
  names(pvals.spearman) = colnames(dat)[-1]
  
  for(mod in names(cc.pearson)) {
    e = pearson.ecdfs[[mod]]
    pvals.pearson[mod] = 1-e(cc.pearson[mod])
    
    e = spearman.ecdfs[[mod]]
    pvals.spearman[mod] = 1-e(cc.spearman[mod])
  }

  mod = "MEbrown"
  fn = file.path(dn.fig, sprintf("pearson-obs-vs-shuffled-hist-%s.pdf",mod))
  pdf(fn, width=5, height=5)
      shuffled = pearson.shuffled[[mod]]
      main = sprintf("%s: Pearson r = %0.3f\n(emp. pval = %0.3f based on %d shuffles)",
                     mod, cc.pearson[mod], pvals.pearson[mod], N)
      hist(shuffled, main=main, xlab="Pearson correlation", breaks=seq(-1,1,0.05), xlim=c(-1,1))
      q95 = quantile(shuffled,0.95)
      lines(c(q95,q95), c(-0.1,1e5), lty="dashed", col="blue", lwd=0.1)
      points(cc.pearson[mod], 0, pch="x", cex=1, col="red")
      try(legend("topright", col=c("blue","red"), pch=c(".", "x"), lty=c("dashed","solid"), lwd=c(0.1,0), c("p=0.05", "observed")))
  dev.off()
  
  df.cor = data.frame(Module = names(cc.pearson), Pearson.cor = cc.pearson, Pearson.empirical.pval = pvals.pearson,
            Spearman.rho = cc.spearman, Spearman.empirical.pval = pvals.spearman)

}


# Table of correlation data ---

df.cor = df.cor %>% 
  mutate(Spearman.FDR = p.adjust(Spearman.empirical.pval, method = "BH")) %>% 
  mutate(Pearson.FDR = p.adjust(Pearson.empirical.pval, method = "BH"))
fwrite(df.cor, file.path(dn.out, "pval-tables_FDR.txt"), sep="\t")
