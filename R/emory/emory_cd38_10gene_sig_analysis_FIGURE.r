library(pROC)

# load gene signature  ---------------------------------------------------------
gene.sig = toupper(fread(file.path(PROJECT_DIR, "generated_data/signatures/CD38_ge_sig.txt"), header=F)[[1]])

box.df = data.frame()
roc.df = data.frame()
df.text = data.frame()

for (year in as.character(c(2008:2011)) ) {
  if(year=="2008") {
    study = "GSE29619"
    load(sprintf(file.path(PROJECT_DIR, "data/Emory/%s.RData"), study,study))
    si = info$year == year & !is.na(info$Subject.ID) & info$time.point=="D0"
  } else {
    study = "GSE74817"
    load(sprintf(file.path(PROJECT_DIR, "data/Emory/%s_healthy.RData"), study,study))
    si = info$year == year & !is.na(info$Subject.ID) & info$time.point=="D0" & info$AGE<=60
  }

# load data
ihl = info$fc.adj.d30[si] %in% c(0,2)
isig = toupper(rownames(dat)) %in% gene.sig
cat(year, sum(isig), "genes\n")
Y = info$fc.adj.d30[si][ihl]
X = dat[isig,si]
if(sum(isig)>1) {
  X = t(scale(t(X)))
  X = colMeans(X)[ihl]
}
w.pv = wilcox.test(X[Y==2], X[Y==0], exact=F, alternative="greater")$p.value

r = roc(Y,X, direction="<", quiet = T)
N_perm = 1000
set.seed(123)
r.null = vector("numeric", N_perm)
for(i in 1:N_perm) {
  r.null[i] = roc(Y, sample(X), direction="<", quiet = T)$auc
}
r.p = sum(r.null > r$auc) / N_perm

r.df = data.frame(Specificity = r$specificities, Sensitivity=r$sensitivities, auc=as.numeric(r$auc), Year=year) %>%
  arrange(Sensitivity, Specificity)
roc.df = rbind(roc.df, r.df)

df.text = rbind(df.text, 
                data.frame(Study = study, Year=year, label.wpv = sprintf("p = %.2f", w.pv),
                           label.auc=sprintf("AUC = %.2f\np = %.2g", r$auc, r.p)))
b.df = data.frame(score=X, Response=Y, Year=year)
box.df = rbind(box.df, b.df)
}


df = box.df %>%
  mutate(Response = factor(ifelse(Response==0, "low","high"), levels=c("low","high")))

source(file.path(PROJECT_DIR, "R/functions/gg_color_hue.r"))
cm = gg_color_hue(4)
p1 = ggplot(df, aes(Response, score, fill=Year)) + geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
  scale_fill_manual(values=c(cm, "black","white")) +
  facet_wrap(~Year, nrow=1) +
  geom_text(data=df.text, aes(label=label.wpv), x=-Inf, y=Inf, col="black", hjust = -0.3, vjust = 1.5) +
  ylab("Relative TGSig score") +
  theme_bw() +
  theme(legend.position="none", 
        panel.border = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
        axis.ticks = element_blank())
p2 = list()
for(year in levels(df.text$Year)) {
  p2[[year]] = ggplot(roc.df %>% dplyr::filter(Year==year), aes(x=Specificity, y=Sensitivity)) + 
    geom_line(size=1) +
    geom_abline(intercept = 1, slope = 1, lty=2) +
    scale_x_reverse() + 
    geom_text(data=df.text %>% dplyr::filter(Year==year), aes(label=label.auc), x=Inf, y=-Inf, hjust=1.1, vjust=-0.1) +
    coord_fixed() + 
    theme_bw() + theme(panel.grid = element_blank())
}

bottom_row <- do.call("plot_grid", c(p2, align = "h", nrow=1))
plot_grid(p1, bottom_row, ncol=1, rel_heights = c(2, 1))

fn.fig = file.path(PROJECT_DIR, "figure_generation/Emory_CD38.10gene.sig_vs_response_4yrs")
ggsave(paste0(fn.fig, ".png"), w=7, h=6)
ggsave(paste0(fn.fig, ".pdf"), w=7, h=6, useDingbats=F)

