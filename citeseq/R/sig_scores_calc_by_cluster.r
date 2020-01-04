library(Seurat)

source("R/functions/get_score.r")

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

df.subj = h1@meta.data %>% mutate(response = str_remove(adjmfc.time, "d0 ")) %>% 
  dplyr::select(subject=sampleid, response) %>% 
  mutate(subject = factor(subject)) %>% 
  distinct() %>% 
  arrange(subject)

dn.out = "results/sig_scores"
dir.create(dn.out, showWarnings = F, recursive = T)

sig.list = readRDS("sig/sig.list.RDS")

for(clustering in 0:3) {
  cat("\n", clustering)
  h1 = SetAllIdent(h1, id = glue::glue("K{clustering}"))
  
  clust.names = levels(h1@ident)

  for(s in names(sig.list)) {
    cat("\n", s, " ")
    sig = sig.list[[ s ]]
    
    score.mat = matrix(nrow = nrow(df.subj), ncol = length(clust.names))
    colnames(score.mat) = clust.names
    rownames(score.mat) = df.subj$subject
    
    for(cl in clust.names) {
      cat(cl, " ")
    
      tobj = SubsetData(h1, ident.use = cl)
      
      dat = tobj@data
      meta = tobj@meta.data %>% 
        dplyr::rename(subject = sampleid) %>% 
        mutate(subject = factor(subject, levels=df.subj$subject))
      
      gi = rownames(dat) %in% sig
      
      dat2 = aggregate(t(as.matrix(dat[gi,,drop=F])), list(subject = meta$subject), mean, drop=F)
      dat2[is.na(dat2)] = 0
      mat = dat2 %>% 
        tibble::column_to_rownames("subject") %>% 
        data.matrix()
    
      score.mat[,cl] = get_score(t(mat))
    } # cl
    # colnames(score.mat) = clust.labels
    score.mat %>% as.data.frame() %>% tibble::rownames_to_column("subject") %>% 
    fwrite(glue::glue("{dn.out}/scores_{s}_{clustering}.txt"), sep = "\t")
  } # sig
} # clustering


