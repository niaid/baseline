library(Seurat)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

sig.list = readRDS("sig/sig.list.RDS")

df.subj = h1@meta.data %>% mutate(response = str_remove(adjmfc.time, "d0 ")) %>% 
  dplyr::select(subject=sampleid, response) %>% 
  mutate(subject = as.numeric(subject)) %>% 
  distinct() %>% 
  arrange(subject)

sigs = c(rep("CD40.act",3), rep("SLE.sig",4), rep("IFN26",2), "TGSig", "LI.M165")
cl.name = c("C3.1.0", "C1", "C6", "C1", "C6", "C7", "C2", "C9", "C2.1.0", "C9", "C9")
cl.level = str_count(cl.name, "\\.")+1
sig.label = glue::glue("{sigs} ({cl.name})")
sig.label = sig.label %>% str_replace("SLE.sig", "SLE-Sig") %>% 
  str_replace("CD40.act", "CD40act") %>% str_replace("IFN26", "IFN-I-DCact")

dn = "results/sig_scores"

df.scores = df.subj
for(s in seq_along(sig.label)) {
  cat(sig.label[s],"\n")
  fn = glue::glue("{dn}/scores_{sigs[s]}_{cl.level[s]}.txt")
  df = fread(fn) %>% dplyr::select("subject",cl.name[s])
  df.scores = left_join(df.scores, df, by="subject")
}
names(df.scores)[-(1:2)] = sig.label


# add CD38++ frequencies in flow
fn.flow = file.path("../generated_data/CHI/CHI_sample_info_2_CD38hi.txt")
df.flow = fread(fn.flow) %>%
  dplyr::filter(time==0) %>%
  dplyr::select(subject, `CD38hi (flow)` = CD38hi)

# add microarray scores
fn.ma = file.path("../generated_data/CHI/CHI_cd38_ge_sig_score_day0.txt")
df.ma = fread(fn.ma) %>% 
  dplyr::select(subject, `TGSig (microarray)` = score)

fn.ma2 = file.path("results/sig_scores/scores_SLE.sig_MA.txt")
df.ma2 = fread(fn.ma2) %>% 
  dplyr::select(subject, `SLE-Sig (microarray)` = SLE.sig)

df.scores = df.scores %>% 
  left_join(df.flow, by="subject") %>% 
  left_join(df.ma, by="subject") %>% 
  left_join(df.ma2, by="subject")
  
fwrite(df.scores, "results/Score_table.txt", sep="\t")

