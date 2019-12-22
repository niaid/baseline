library(effects)

source(file.path(PROJECT_DIR, "R/brown-mod-minus-leading-edge-combined/input.R"))

dn.fig = file.path(PROJECT_DIR, "figure_generation/SLE-Sig")
dir.create(dn.fig, showWarnings = F)


flu[, Study:=as.factor(Study)]

##
## one predictor: only brown module minus leading edge genes
##
fit1 = glm(HighResponse ~ brown.minus.LE + Study, data=flu, family="binomial")
summary(fit1)

fn = file.path(dn.fig, "4flu-brown-MINUS-LE-effects-one-predictor-bars.pdf")
pdf(fn, width=5, height=5)
  ylab = "Probability of high responder"
  main = sprintf("p (brown_minus_LE)=%0.3f",
                 coef(summary(fit1))[2,4])
  eff = Effect("brown.minus.LE", fit1, xlevels=list(brown.minus.LE=seq(-0.5, 0.5, 0.25)))
  barplot(summary(eff)$effect, ylim=c(0,1), xlab="Brown module excluding LE genes", ylab=ylab, width=1, space=0, las=1,
          main=main)
  with(summary(eff),
       gplots::plotCI(seq(0.5, 4.5, 1), effect, uiw = upper - effect, liw = effect - lower, add=T, pch=NA,
                      gap=0, type="p")
       )
dev.off()

