library(pROC)
source(file.path(PROJECT_DIR, "R/functions/load_sig.r"))
source(file.path(PROJECT_DIR, "R/functions/get_score.r"))

# load gene expression data
fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix_gene.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info_2.txt")
info = fread(fn.si)

ng=10
fn.cd38.cor = file.path(PROJECT_DIR, "generated_data", "CHI", "robust_corr_genes.txt")
gene.sig = load_sig(fn.cd38.cor, "cor.mean.sd.ratio", ntop=ng)

fn.rand = file.path(PROJECT_DIR, "generated_data/CHI/gene_sig_random_500_from_all.txt")
if(!file.exists(fn.rand)) {
  fn.rand = file.path(PROJECT_DIR, "data/random_sig/gene_sig_random_500_from_all.txt")
}
if(!file.exists(fn.rand)) {
  stop("No file with random signatures found in data or generated_data\nRun ./R/chi_signature_analysis/rand.sig_generate.r")
}

r.sig.df = read.table(fn.rand, sep="\t",header=F, row.names = 1)

roc.df = data.frame()
rand.df = data.frame()
N_rand = 500
N_rand = min(N_rand, ncol(r.sig.df))

for (time.pt in c(0,-7,70)) {
  cat(time.pt,"")
  
  day.pt = paste0("day ", time.pt)
  # test signature
  gi = toupper(rownames(dat)) %in% toupper(gene.sig)
  sum(gi)

  iday = info$time==time.pt
  ihl = info$adjMFC_class[iday] %in% c(0,2)
  X = get_score(dat[gi,iday])[ihl]
  Y = info$adjMFC_class[iday][ihl]
  
  r = roc(Y,X, direction="<", quiet = T)
  r.df = data.frame(day=day.pt, AUC = as.numeric(r$auc))
  roc.df = rbind(roc.df, r.df)
  
  # test AUC in random signatures
  for (i in 1:N_rand) {
    r.sig = load_sig(r.sig.df, i, ntop=ng)
    gi = toupper(rownames(dat)) %in% r.sig
    sum(gi)
    
    X = get_score(dat[gi,iday])[ihl]

    r = roc(Y,X, direction="<", quiet = T)
    r.df = data.frame(day=day.pt, AUC = as.numeric(r$auc))
    rand.df = rbind(rand.df, r.df)
  }
  
}

fn.rand.rds = file.path(PROJECT_DIR, "generated_data/CHI/chi_random.sig.500_3times.rds")
saveRDS(rand.df, "chi_random.sig.500_3times.rds")

# fn.rand.rds = file.path(PROJECT_DIR, "generated_data/CHI/chi_random.sig.500_3times.rds")
# rand.df = readRDS(fn.rand.rds)

roc.df = roc.df %>% 
  mutate(day = case_when(
    day == "day 0" ~ "Baseline 1 (day 0)",
    day == "day -7" ~ "Baseline 2 (day -7)",
    day == "day 70" ~ "Baseline 3 (day 70)"
  ) %>% factor())
rand.df = rand.df %>% 
  mutate(day = case_when(
    day == "day 0" ~ "Baseline 1 (day 0)",
    day == "day -7" ~ "Baseline 2 (day -7)",
    day == "day 70" ~ "Baseline 3 (day 70)"
  ) %>% factor())

rand.df.p = rand.df %>% 
  rename(rAUC = AUC) %>% 
  left_join(roc.df, by="day") %>% 
  group_by(day, AUC) %>% 
  summarise(p=sum(rAUC>AUC)/nrow(.)) %>% 
  ungroup()

ggplot(rand.df, aes(AUC, fill=day)) + 
  geom_density(col=NA, alpha=1) + 
  facet_wrap(~day,ncol=1) +
  scale_fill_manual(values=c("grey65", "grey80","grey50")) +
  xlim(c(0,1)) +
  geom_vline(data=roc.df, aes(xintercept=AUC), col="red", lty=2) +
  geom_text(data=rand.df.p, aes(label=sprintf("p = %.1g",p)), x=-Inf, y=Inf, hjust=-0.1, vjust = 1.1, size=4) +
  theme_bw() + guides(fill="none") + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=12))

fn.fig = file.path(PROJECT_DIR, "figure_generation", 
                   "CHI_CD38_10genes_random_sig")
ggsave(paste0(fn.fig,".png"), w=4,h=7)
ggsave(paste0(fn.fig,".pdf"), w=4,h=7)
