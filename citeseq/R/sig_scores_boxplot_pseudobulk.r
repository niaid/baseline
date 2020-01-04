library(ggsignif)

df.test = fread("results/test_sig.genes_in_clusters_high_vs_low_responders.txt") %>% 
  dplyr::filter(w.pv < 0.05, test.dir == "IN", test == "SigScore", 
                clustering == "K0") %>% 
  mutate(sig.name = sig %>% 
          str_replace("SLE.sig", "SLE-Sig") %>% 
          str_replace("CD40.act", "CD40act") %>% 
          str_replace("IFN26", "IFN-I-DCact")) %>% 
  mutate(sig = factor(sig)) %>% 
  arrange(sig, cluster)

clustering = 0
sigs = levels(df.test$sig)
sig.label = sigs %>% 
  str_replace("SLE.sig", "SLE-Sig") %>% 
  str_replace("CD40.act", "CD40act") %>% 
  str_replace("IFN26", "IFN-I-DCact") %>% 
  setNames(sigs)


dn.fig = "figures/sig_test_boxplots"
dir.create(dn.fig, showWarnings = F, recursive = T)

for(s in sigs) {
  cat(s, "\n")
fn.score = glue::glue("results/sig_scores/scores_{s}_{clustering}.txt")
df.sc = fread(fn.score, data.table = F) %>% 
  mutate(subject = factor(subject)) %>% 
  gather("cluster","score", -subject) %>% 
  left_join(df.subj, by="subject") %>% 
  # mutate(clustern = 9-as.numeric(str_remove(cluster,"C"))) %>%
  mutate(response = factor(response, levels=c("low","high"))) %>%
  arrange(cluster, response) %>% 
  mutate(resp_cl = ifelse(response=="low", -0.15, 0.15))

df.txt = df.test %>% 
  dplyr::filter(sig == s) %>% 
  mutate(x=-Inf, y=Inf, label=glue::glue("p = {format(w.pv, digits=2)}"))
yrng = c(-1, 1.55)
cm = pals::glasbey(13)[c(1:3,5:9,13,11)]
dotsz = diff(yrng)/(max(df.sc$score) - min(df.sc$score))
p = ggplot(df.sc, aes(response, score, group=response)) +
  coord_flip(ylim=yrng) + 
  geom_boxplot(width=0.6, aes(fill=cluster), position = position_dodge2(preserve = "total"), lwd = 0.3, alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_dotplot(aes(fill = response), binaxis = "y", stackdir = "center", dotsize = dotsz, 
               position=position_dodge(), show.legend = F) +
  scale_fill_manual(values = c("black", "white", "white")) +
  geom_text(aes(x, y, label=label), col="red", hjust=1, vjust=-0.5, data=df.txt, inherit.aes = F) +
  xlab("") + ylab(paste(sig.label[s], "score")) +
  theme_classic()

ggsave(glue::glue("{dn.fig}/sig_scores_box_{s}_K{clustering}.png"), plot=p, w=2.5, h=1.2)
ggsave(glue::glue("{dn.fig}/sig_scores_box_{s}_K{clustering}.pdf"), plot=p, w=2.5, h=1.2, useDingbats=F)
}
