library(pROC)
library(cowplot)
source(file.path(PROJECT_DIR, "R/functions/get_score.r"))
source(file.path(PROJECT_DIR, "R/functions/load_sig.r"))

# load CD38 signature genes
fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_genes.txt")
gene.sig = load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=10)

# load gene expression data
fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix_gene.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info_2_CD38hi.txt")
info = fread(fn.si, data.table=F) %>% 
  mutate(subject = as.character(subject)) %>% 
  mutate(Response = factor(Response, levels=c("low","high")))

roc.df = data.frame()
df = data.frame()
df.text = data.frame()
N_perm = 1000

for (tp in c(0,-7,70)) {
  si = with(info, time == tp)# & !is.na(CD38hi))
  sum(si)
  tdat = dat[,si]
  tinfo = info[si,]
  
  day.pt = case_when(
    tp == 0 ~ "Baseline 1 (day 0)",
    tp == -7 ~ "Baseline 2 (day -7)",
    tp == 70 ~ "Baseline 3 (day 70)"
  )
  
  gi = toupper(rownames(tdat)) %in% toupper(gene.sig)
  sum(gi)
  
  ihl = tinfo$Response %in% c("low","high")
  X = get_score(tdat[gi,])[ihl]
  Y = tinfo$Response[ihl]
  
  r = roc(Y,X, direction="<", quiet = T)
  set.seed(123)
  r.null = vector("numeric", N_perm)
  for(i in 1:N_perm) {
    r.null[i] = roc(Y, sample(X), direction="<", quiet = T)$auc
  }
  r.p = sum(r.null > r$auc) / N_perm
  r.df = data.frame(day.label=day.pt, Specificity = r$specificities, Sensitivity=r$sensitivities) %>%
    arrange(Sensitivity)
  
  roc.df = rbind(roc.df, r.df)
  
  w.pv = wilcox.test(X[Y=="high"], X[Y=="low"], exact=F, alternative = "greater")$p.value
  
  df.text = rbind(df.text, 
                  data.frame(day.label = day.pt, label.auc=sprintf("AUC = %.2f\np = %.2g", r$auc, r.p), 
                             label.w.pv=sprintf("p = %.2g",w.pv)))
  
  tmp = data.frame(day = paste0("day ", tp), day.label=day.pt, subject=tinfo$subject[ihl], score=X, Response=factor(Y))
  df = rbind(df, tmp)
}

df.text$x = rep(0.25, 3)
df.text$y = rep(0.1, 3)

roc.df = roc.df %>% 
  mutate(day = factor(day.label))
df = df %>% 
  mutate(day.label = factor(day.label))

# output scores
fn.out = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_cd38_ge_sig_score.txt")
fwrite(df %>% dplyr::select(-day.label), fn.out, sep="\t", quote=T)

fn.out = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_cd38_ge_sig_score_day0.txt")
fwrite(df %>% dplyr::select(-day.label) %>% dplyr::filter(day=="day 0"), fn.out, sep="\t", quote=T)

# Generate boxplots
p1 = ggplot(df, aes(x=Response, y=score, group=Response)) +
  geom_boxplot(aes(fill=day.label), alpha=1, outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
  facet_wrap(~day.label, nrow=1) +
  scale_fill_manual(values=c("grey70", "grey90","grey50", "black", "white"), name="Response",
                    breaks=c(0,2), labels=c("low","high")) +
  geom_text(data = df.text, aes(label=label.w.pv), x=1.5, y=Inf,vjust=1.1, hjust=0.5, size=4, inherit.aes = F) +
  xlab("Response") + ylab("Baseline signature score") +
  coord_cartesian(xlim=c(0.8, 2.2)) +
  theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
        axis.ticks = element_blank())

p2 = list()
for(tm in levels(df.text$day.label)) {
  p2[[tm]] = ggplot(roc.df %>% dplyr::filter(day.label==tm), aes(x=Specificity, y=Sensitivity)) + 
    geom_line(size=1) +
    geom_abline(intercept = 1, slope = 1, lty=2) +
    scale_x_reverse() + 
    geom_text(data=df.text %>% dplyr::filter(day.label==tm), aes(x=x,y=y,label=label.auc)) +
    coord_fixed() + 
    theme_bw() + theme(panel.grid = element_blank())
}

bottom_row <- do.call("plot_grid", c(p2, align = "h", nrow=1))
plot_grid(p1, bottom_row, ncol=1, rel_heights = c(2, 1))

fn.fig = file.path(PROJECT_DIR, "figure_generation", "CD38.10gene.sig_vs_response_3time_p")
ggsave(paste0(fn.fig, ".png"), w=7,h=7)
ggsave(paste0(fn.fig, ".pdf"), w=7,h=7)

