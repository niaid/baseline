library(gplots)
library(car)
library(pROC)
library(data.table)
library(MASS)
library(effects)


source(file.path(PROJECT_DIR, "R/recalc-and-check-10gene-CD38plusSig/input.R"))

dn.fig = file.path(PROJECT_DIR, "figure_generation/SLE-Sig")
dir.create(dn.fig, showWarnings = F)

summary(st)

partial.residuals=F
###########
if(1) {
  ## TWO classes:
  st2 = st[which(Response != "1")]
  st2[, Resp:=as.factor(Response)]
  st2[, Study:=as.factor(Study)]

  ##
  ## predictor: using CD38 score only
  ##
  fit1 = glm(Resp ~ CD38_score + Study, family="binomial", data=st2)
  summary(fit1)
  fit1.coef=coef(summary(fit1))
  fit1.coef

  fn = file.path(dn.fig, "4flu-logistic-regression-effects-CD38sig-bars.pdf")
  pdf(fn, width=5, height=5)
    p = summary(fit1)$coef[2,4]
    main=sprintf("logistic regression (pval=%0.3g)", p)
         ylab = "Probability of high responder"
         eff = Effect("CD38_score", fit1, xlevels=list(CD38_score=seq(-1,1,0.5)))
    with(summary(eff), {
         barplot(effect, main=main,  width=1, space=0, ylim=c(0,1), las=1,
                 xlab="TGSig score", ylab=ylab)
         plotCI(seq(0.5,4.5,1), effect, uiw=upper-effect, liw=effect-lower, add=T, gap=0, type="p", pch=NA)
                  })
  dev.off()
  
  ##
  ## using brown module only
  ##
  fit1 = glm(Resp ~ brown.L.E.eigengene + Study, family="binomial", data=st2)
  summary(fit1)
  fit1.coef=coef(summary(fit1))
  fit1.coef
  
  fn = file.path(dn.fig, "4flu-logistic-regression-effects-brown-leading-edge-eigengene-bars.pdf")
  pdf(fn, width=5, height=5)
    p = summary(fit1)$coef[2,4]
    main=sprintf("logistic regression (pval=%0.3g)", p)
    ylab = "Probability of high responder"
    eff = Effect("brown.L.E.eigengene", fit1, xlevels=list(brown.L.E.eigengene=seq(-1,1,0.5)))
    with(summary(eff), {
         barplot(effect, main=main,  width=1, space=0, ylim=c(0,1), las=1,
                 xlab="SLE-Sig score", ylab=ylab)
         plotCI(seq(0.5,4.5,1), effect, uiw=upper-effect, liw=effect-lower, add=T, gap=0, type="p", pch=NA)
                  })
  dev.off()

  
  ##
  ## using both CD38 score and brown module 
  ##
  fit2 = glm(Resp ~ CD38_score + brown.L.E.eigengene + Study, family="binomial",
              data=st2)
  summary(fit2)
  fit2.coef=coef(summary(fit2))
  fit2.coef
  fn = file.path(dn.fig, "4flu-logistic-regression-effects-two-predictors-showing-CD38sig-bars.pdf")
  pdf(fn, width=5, height=5)
    eff = Effect("CD38_score", fit2, xlevels=list(CD38_score=seq(-1,1,0.5)))
    main = sprintf("pval (TGSig) = %0.5f", summary(fit2)$coef[2,4])
    barplot(summary(eff)$effect, ylim=c(0,1), xlab="TGSig score", ylab=ylab, width=1, space=0, las=1,
            main=main)
    with(summary(eff), 
         gplots::plotCI(seq(0.5, 4.5, 1), effect, uiw = upper - effect, liw = effect - lower, add=T, pch=NA,
                        gap=0, type="p")
         )
  dev.off()
  
  fn = file.path(dn.fig, "4flu-logistic-regression-effects-two-predictors-showing-brown-leading-edge-eigengene-bars.pdf")
  pdf(fn, width=5, height=5)
    eff = Effect("brown.L.E.eigengene", fit2, xlevels=list(brown.L.E.eigengene=seq(-1,1,0.5)))
    main = sprintf("pval (SLE-Sig) = %0.5f", summary(fit2)$coef[3,4])
    barplot(summary(eff)$effect, ylim=c(0,1), xlab="SLE-Sig score", ylab=ylab, width=1, space=0, las=1,
            main=main)
    with(summary(eff), 
         gplots::plotCI(seq(0.5, 4.5, 1), effect, uiw = upper - effect, liw = effect - lower, add=T, pch=NA,
                        gap=0, type="p")
         )
  dev.off()
  
}

