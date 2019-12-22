source(file.path(PROJECT_DIR, "R/functions/get_score.r"))
source(file.path(PROJECT_DIR, "R/functions/load_sig.r"))
source(file.path(PROJECT_DIR, "R/functions/gg_color_hue.r"))

# load CD38 signature genes
fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_genes.txt")
gene.sig = load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=10)

# load gene expression data
fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix_gene.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info_2_CD38hi.txt")
info = fread(fn.si, data.table=F)
info = info %>% mutate(Response = factor(Response, levels=c("low","middle","high")))

df = data.frame()
df.text = data.frame()

for (tp in c(0,-7,70)) {
  si = with(info, time == tp & !is.na(Response))
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
  
  X = get_score(tdat[gi,])
  Y = tinfo$gender
  
  w.pv = wilcox.test(X[Y=="F"], X[Y=="M"], exact=F)$p.value
  df.text = rbind(df.text, 
                  data.frame(day.label = day.pt, 
                             label.w.pv=sprintf("p = %.2g",w.pv)))
  
  tmp = data.frame(day = paste0("day ", tp), day.label=day.pt, subject=tinfo$subject, 
                   score=X, Gender=factor(Y), Response=factor(tinfo$Response))
  df = rbind(df, tmp)
}

df.text$x = rep(0.25, 3)
df.text$y = rep(0.1, 3)

df = df %>% 
  mutate(day.label = factor(day.label))

# Generate boxplots
clr = gg_color_hue(2)
ggplot(df, aes(x=Gender, y=score)) +
  geom_boxplot(aes(fill=Gender, alpha=day.label), outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
  facet_wrap(~day.label, nrow=1) +
  scale_fill_manual(values=c(clr[1], "black","white", clr[2], "grey50"), na.value="grey80") +
  scale_alpha_manual(values = c(0.5,0.2,0.8)) +
  geom_text(data = df.text, aes(label=label.w.pv), x=1.5, y=Inf,vjust=1.1, hjust=0.5, size=4, inherit.aes = F) +
  xlab("Response") + ylab("TGSig score") +
  coord_cartesian(xlim=c(0.8, 2.2)) +
  theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
        axis.ticks = element_blank())

fn.fig = file.path(PROJECT_DIR, "figure_generation", glue::glue("CHI_TGSig_vs_gender_BL_all_subjects"))
ggsave(paste0(fn.fig, ".png"), w=7,h=4)
ggsave(paste0(fn.fig, ".pdf"), w=7,h=4, useDingbats=F)

