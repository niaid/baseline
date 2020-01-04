library(Seurat)

source("R/functions/TestGeneSig.r")

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

# load signatures
sig.list = readRDS("sig/sig.list.RDS")

df.test = fread("results/test_sig.genes_in_clusters_high_vs_low_responders.txt") %>% 
  dplyr::filter(w.pv < 0.05, test.dir == "IN", test == "SigScore", cluster != "ALL") %>% 
  mutate(sig_cluster = paste(sig, cluster, sep="_")) %>% 
  mutate(sig = factor(sig, levels = names(sig.list))) %>% 
  arrange(clustering, sig, cluster)

df.subj = h1@meta.data %>% mutate(response = str_remove(adjmfc.time, "d0 ")) %>% 
  dplyr::select(subject=sampleid, response) %>% 
  mutate(subject = factor(subject)) %>% 
  distinct() %>% 
  arrange(subject)


df.out = data.frame()

clusters.out = c(as.list(paste0("C",0:9)), 
                      list(c("C1","C6"), c("C1","C6","C9","C2.1.0"), c("C3.1.0"), c("C9","C2.1.0")))
df.cells = h1@meta.data %>% dplyr::select(K1,K3) %>% 
  tibble::rownames_to_column("cell")

for(cl.out in clusters.out) {
  cat(cl.out, "\n")
  tobj = SubsetData(h1, cells.use = df.cells %>% dplyr::filter(!(K1 %in% cl.out | K3 %in% cl.out)) %>% pull("cell"))
  cl.out.name = paste(cl.out, collapse = ", ")
  for(s in names(sig.list)[c(1,2,6)]) {
      df.out = rbind(df.out,
                      TestGeneSig(tobj, sig.list[s], "sampleid", df.subj, 
                                  "response", c("low","high"), test.to.print = NA) %>% 
                        mutate(cluster = "Filtered", test.dir="OUT", label = cl.out.name)
      )
  } # sig
  cat("\n")
} # clusters dropped

fwrite(df.out, "results/test_sig.genes_in_clusters_high_vs_low_responders_OUT.txt", sep="\t")


sig.in = c("TGSig", "SLE.sig", "CD40.act")
res = fread("results/test_sig.genes_in_clusters_high_vs_low_responders.txt") %>% 
  dplyr::filter(cluster=="pseudo-bulk", test=="SigScore") %>% 
  dplyr::select(-clustering) %>% 
  mutate(label="")
res.out = fread("results/test_sig.genes_in_clusters_high_vs_low_responders_OUT.txt") 
  
res = rbind(res, res.out) %>% 
  dplyr::filter(sig %in% sig.in, test=="SigScore")

df.test.fig = res %>% 
  mutate(sig = ifelse(sig=="SLE.sig", "SLE-Sig", ifelse(sig=="CD40.act", "CD40act", sig))) %>%
  mutate(label = ifelse(label=="", "pseudo-bulk", label)) %>% 
  mutate(cluster = fct_inorder(as.character(cluster))) %>% 
  mutate(label = fct_inorder(label) %>% fct_rev()) %>% 
  mutate(sig = fct_inorder(as.character(sig))) %>% 
  mutate(signif = case_when(
    w.pv<=0.001 ~ "***",
    w.pv<=0.01 ~ "**",
    w.pv<=0.05 ~ "*",
    TRUE ~ ""
  ))
ggplot(df.test.fig, aes(sig, label, fill=-log10(w.pv))) +
  geom_tile() +
  geom_text(aes(label=signif), size=6) +
  scale_fill_gradient(low="white", high="red", name="-log10(p)") +
  xlab(NULL) + 
  ylab("Clusters dropped out") +
  theme_minimal()

dn.fig = "figures/drop_out_test"
dir.create(dn.fig, showWarnings = F, recursive = T)
fn.fig = file.path(dn.fig, "clusters_drop_out_heatmap")
ggsave(paste0(fn.fig, ".png"), w=3.9, h=3.5)
ggsave(paste0(fn.fig, ".pdf"), w=3.9, h=3.5, useDingbats=F)
