library(ggsignif)

df.test = fread("results/test_sig.genes_in_clusters_high_vs_low_responders.txt") %>% 
  dplyr::filter(w.pv < 0.05, test.dir == "IN", test == "SigScore", 
                clustering == "K1") %>% 
  dplyr::mutate(sig.name = sig %>% 
          str_replace("SLE.sig", "SLE-Sig") %>% 
          str_replace("CD40.act", "CD40act") %>% 
          str_replace("IFN26", "IFN-I-DCact")) %>% 
  dplyr::mutate(clustern = 9-as.numeric(str_remove(cluster,"C"))) %>%
  dplyr::mutate(sig = factor(sig)) %>% 
  dplyr::arrange(sig, cluster)

clustering = 1
sigs = levels(df.test$sig)
sig.label = sigs %>% str_replace("SLE.sig", "SLE-Sig") %>% 
  str_replace("CD40.act", "CD40act") %>% str_replace("IFN26", "IFN-I-DCact") %>% 
  setNames(sigs)


dn.fig = "figures/sig_test_boxplots"
dir.create(dn.fig, showWarnings = F, recursive = T)

for(s in sigs) {
  cat(s, "\n")
fn.score = glue::glue("results/sig_scores/scores_{s}_{clustering}.txt")
df.sc = fread(fn.score, data.table = F) %>% 
  dplyr::mutate(subject = factor(subject)) %>% 
  gather("cluster","score", -subject) %>% 
  left_join(df.subj, by="subject") %>% 
  dplyr::mutate(clustern = 9-as.numeric(str_remove(cluster,"C"))) %>%
  dplyr::mutate(response = factor(response, levels=c("low","high"))) %>%
  dplyr::arrange(cluster, response) %>% 
  dplyr::mutate(resp_cl = ifelse(response=="low", clustern-0.15, clustern+0.15))

df.txt = df.test %>% 
  dplyr::filter(sig == s) %>% 
  dplyr::mutate(x = clustern, y = 1, label = "*")
yrng = c(-1, 1.55)
cm = pals::glasbey(13)[c(1:3,5:9,13,11)]
dotsz = diff(yrng)/(max(df.sc$score) - min(df.sc$score))
p = ggplot(df.sc, aes(resp_cl, score, group=resp_cl)) +
  coord_flip(ylim=yrng) + 
  geom_boxplot(width=0.6, aes(fill=cluster), position = position_dodge2(preserve = "total"), lwd = 0.3, alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_dotplot(aes(fill = response), binaxis = "y", stackdir = "center", dotsize = dotsz, 
               position=position_dodge(), show.legend = F) +
  scale_fill_manual(values = c(cm,"black", "white")) +
  scale_x_continuous(breaks=0:10 ,labels = c(paste0("C",9:0) ,"All")) +
  xlab("Cell cluster") + ylab(paste(sig.label[s], "score")) +
  theme_classic()
if(nrow(df.txt)>0) {
  p = p + ggsignif::geom_signif(y_position = rep(1,nrow(df.txt)), textsize = 8, col="red",
                      xmin = df.txt$x-0.15, xmax = df.txt$x+0.15, 
                      annotation=df.txt$label, tip_length = 0.01)
}

ggsave(glue::glue("{dn.fig}/sig_scores_box_{s}_K{clustering}.png"), plot=p, w=2, h=6)
ggsave(glue::glue("{dn.fig}/sig_scores_box_{s}_K{clustering}.pdf"), plot=p, w=2, h=6, useDingbats=F)
}
