library(lme4)
source("R/functions/get_score.r")

fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle.txt")
info = fread(fn.si, data.table = F)

info = info %>% 
  mutate(DA=factor(DA, levels=c("low","mid","high")), 
         SUBJECT = factor(SUBJECT, levels=
                            unique(SUBJECT[order(as.numeric(sub("SLE-","",SUBJECT)))])))

info = info %>% 
  mutate(DA=factor(DA, levels=c("low","mid","high")), 
         DAn = as.numeric(DA),
         SUBJECT = factor(SUBJECT, levels=
                            unique(SUBJECT[order(as.numeric(sub("SLE-","",SUBJECT)))])))

info = info %>% 
  mutate(DA2 = ifelse(DAn>1, "high", "low")) %>% 
  group_by(SUBJECT) %>% 
  mutate(DAn.min = min(DAn, na.rm=t), DAn.max = max(DAn, na.rm=t)) %>% 
  ungroup()

si = info$DAn.min==1 & info$DAn.max>1 & (info$DAn==info$DAn.min | info$DAn==info$DAn.max)
sum(si)
subj = info$SUBJECT[si] %>% as.character() %>% unique()


treatment.fields = c("STEROID_IV_CATEGORY","CYCLOPHOSPHAMIDE_CATEGORY",
                     "ORAL_STEROIDS_CATEGORY","MYCOPHENOLATE_CATEGORY",
                     "HYDROXYCHLOROQUINE_CATEGORY")
info = info %>% mutate_at(.vars=treatment.fields, .funs = factor)

fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene_sle.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()

dat = dat[,si]
info = info[si,]

pb.file = "PB_DC.M4.11_ge_sig.txt"
ppb.label = "PB.DC"

fn.pg = file.path(PROJECT_DIR, "data", "SLE", "phenotypes", "SLE_SUBJECT_PG.txt")
df.pg = read_tsv(fn.pg)


fn.sig = file.path(PROJECT_DIR, "generated_data", "signatures", pb.file)
pb.genes = fread(fn.sig, header = F) %>% unlist(use.names=F)

gi = toupper(rownames(dat)) %in% toupper(pb.genes)
sum(gi)
info = mutate(info, PB=get_score(dat[gi,]))

pbm = info %>% 
  group_by(SUBJECT, DA2) %>% 
  summarise(PB = mean(PB, na.rm=T)) %>% 
  ungroup()
pb.fc = pbm %>% 
  spread("DA2","PB") %>% 
  mutate(dPB = high-low)


form = paste0("PB"," ~ 1 + RACE + ",
              paste(treatment.fields,sep="",collapse=" + "), " + (SLEDAI|SUBJECT)")
m.lmer = lmer(as.formula(form), data=info)
mre = ranef(m.lmer, condVar=T)
x = mre$SUBJECT
pv = attr(x, "postVar")
se = unlist(lapply(1:ncol(x), function(i) sqrt(pv[i, i, ])))
mre.df = data.frame(y=unlist(x),
                    se=se,
                    ci=1.96*se,
                    nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                    SUBJECT=factor(rep(rownames(x), ncol(x)), levels=rownames(x)),
                    group=gl(ncol(x), nrow(x), labels=names(x))) %>%
  dplyr::filter(group %in% c("SLEDAI")) %>%
  mutate(z = y/ci)
mre.z = mre.df %>%
  dplyr::select(SUBJECT,group,z) %>%
  spread(group,z) %>% 
  dplyr::rename(DaCP=SLEDAI) 

df.out = left_join(mre.z, df.pg, by="SUBJECT")

df.out = df.out %>% 
  left_join(pb.fc, by="SUBJECT")

PG.group.names = c("SLE patients with plasmablast\nsignature during flare", "Other SLE patients")

df.fig = df.out %>% 
  dplyr::filter(!is.na(PG)) %>% 
  mutate(PG2 = ifelse(PG %in% 2:3, "2/3", ifelse(PG %in% 4, "4", "other"))) %>% 
  mutate(PG.group = ifelse(PG %in% 2:4, PG.group.names[1], PG.group.names[2])) %>% 
  mutate(PG = factor(PG), PG.group = factor(PG.group, levels=PG.group.names))
df.cc = df.fig %>% split(.$PG.group) %>% 
  map(~cor.test(.$dPB, .$DaCP, method="pearson")) %>% 
  map_df(~data.frame(label=sprintf("r = %.3f\np = %.2g",.$estimate,.$p.value),
                     x=-Inf, y=Inf)) %>% 
  mutate(PG.group = levels(df.fig$PG.group) %>% factor())

df.fig2 = df.fig %>% dplyr::filter(PG %in% 2:3)
cor.test(~dPB+DaCP, data = df.fig2, method="pearson") 

ggplot(df.fig %>% dplyr::filter(PG %in% 2:4), aes(dPB, DaCP)) +
  geom_point(aes(fill=PG), size=3, shape=21, col="white", show.legend = F) + 
  scale_fill_manual(values = c("darkblue","darkblue","grey70")) +
  geom_text(data=df.cc[1,], aes(label=label), x=-Inf, y=Inf, vjust=1.5, hjust=-0.1, col="black") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) + theme_bw() + theme(legend.key = element_blank()) + 
  xlab("Delta between mean PB score") + 
  ylab("DaCP") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(10,"mm") )

fn.fig = file.path(PROJECT_DIR, "figure_generation/SLE_dPB_DaCP_correlation")
ggsave(paste0(fn.fig,".png"), w=4,h=3)
ggsave(paste0(fn.fig,".pdf"), w=4,h=3, useDingbats=F)
