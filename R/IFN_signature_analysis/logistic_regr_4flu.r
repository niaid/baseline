fn.st = file.path(PROJECT_DIR, "generated_data/brown-module-leading-genes-low-mid-high-eigengene/st.rds")
st = readRDS(fn.st) %>% 
  dplyr::rename(subject = Sample, CD38 = CD38_score) %>% 
  mutate(Response = factor(Response), Study = factor(Study))


# CHI IFN26
fn.chi = file.path(PROJECT_DIR, "generated_data", "IFN26", "CHI_IFN26_ge_sig_score.txt")
df.chi = fread(fn.chi) %>% 
  dplyr::filter(day == "day 0") %>% 
  dplyr::select(subject, IFN26 = score)

# HIPC IFN26
fn.hipc = file.path(PROJECT_DIR, "generated_data", "IFN26", "HIPC_IFN26_ge_sig_score.txt")
df.hipc = fread(fn.hipc) %>% 
  dplyr::select(subject, IFN26 = score)

# CHI TGSig
fn.cd38 = file.path(PROJECT_DIR, "generated_data/CHI/CHI_cd38_ge_sig_score_day0.txt")
df.cd38 = fread(fn.cd38)

df = rbind(df.chi, df.hipc)

DF = st %>% 
  inner_join(df, by="subject")

dn.fig = file.path(PROJECT_DIR, "figure_generation/IFN26")


### logstic regression - 1 predictor

fit1 = glm(Response ~ IFN26 + Study, family="binomial", data=DF)
summary(fit1)
# fit1.coef=coef(summary(fit1))
# fit1.coef

fn.fig = file.path(dn.fig, "4flu-logistic-regression-1pred_effects-IFN26-bars")
p = summary(fit1)$coef[2,4]
main=sprintf("p = %0.3g", p)
ylab = "Probability of high responder"
eff = Effect("IFN26", fit1, xlevels=list(IFN26=seq(-1,1,0.5)))
with(summary(eff), {
  barplot(effect, main=main,  width=1, space=0, ylim=c(0,1), las=1,
          xlab="IFN-I-DCact score", ylab=ylab)
  plotCI(seq(0.5,4.5,1), effect, uiw=upper-effect, liw=effect-lower, add=T, gap=0, type="p", pch=29)
})
dev.copy(png, paste0(fn.fig, ".png"))
dev.off()
dev.copy(pdf, paste0(fn.fig, ".pdf"))
dev.off()


### logstic regression - 2 predictors

fit1 = glm(Response ~ CD38 + IFN26 + Study, family="binomial", data=DF)
summary(fit1)
# fit1.coef=coef(summary(fit1))
# fit1.coef

fn.fig = file.path(dn.fig, "4flu-logistic-regression-2pred_effects-CD38-bars")
p = summary(fit1)$coef[2,4]
main=sprintf("p = %0.3g", p)
ylab = "Probability of high responder"
eff = Effect("CD38", fit1, xlevels=list(CD38=seq(-1,1,0.5)))
with(summary(eff), {
  barplot(effect, main=main,  width=1, space=0, ylim=c(0,1), las=1,
          xlab="TGSig score", ylab=ylab)
  plotCI(seq(0.5,4.5,1), effect, uiw=upper-effect, liw=effect-lower, add=T, gap=0, type="p", pch=NA)
})
dev.copy(png, paste0(fn.fig, ".png"))
dev.off()
dev.copy(pdf, paste0(fn.fig, ".pdf"))
dev.off()

fn.fig = file.path(dn.fig, "4flu-logistic-regression-2pred_effects-IFN26-bars")
p = summary(fit1)$coef[3,4]
main=sprintf("p = %0.3g", p)
ylab = "Probability of high responder"
eff = Effect("IFN26", fit1, xlevels=list(IFN26=seq(-1,1,0.5)))
with(summary(eff), {
  barplot(effect, main=main,  width=1, space=0, ylim=c(0,1), las=1,
          xlab="IFN-I-DCact score", ylab=ylab)
  plotCI(seq(0.5,4.5,1), effect, uiw=upper-effect, liw=effect-lower, add=T, gap=0, type="p", pch=NA)
})
dev.copy(png, paste0(fn.fig, ".png"))
dev.off()
dev.copy(pdf, paste0(fn.fig, ".pdf"))
dev.off()

