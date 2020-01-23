library(Seurat)
source("R/functions/TestGeneSig.r")

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

dn.out = "results"
dir.create(dn.out, showWarnings = F, recursive = T)

df.subj = h1@meta.data %>% dplyr::mutate(response = str_remove(adjmfc.time, "d0 ")) %>% 
  dplyr::select(subject=sampleid, response) %>% 
  dplyr::mutate(subject = factor(subject)) %>% 
  distinct() %>% 
  dplyr::arrange(subject)

# load signatures
sig.list = readRDS("sig/sig.list.RDS")

df.test = data.frame()

for(clustering in paste0("K",0:3)) {
  cat(clustering, "\n")
  h1 = SetAllIdent(h1, id = clustering)
  if(clustering != "K0") {
    cluster.names = h1@meta.data %>% pull(clustering) %>% unique() %>% sort()
  } else {
    cluster.names = "pseudo-bulk"
  }
  for(test.dir in c("IN")) {
    cat(test.dir, "\n")
    for (ccl.name in cluster.names) {
      if(clustering == "K0") {
        if(test.dir != "IN") next()
        tobj = h1
      } else {
        if(test.dir == "IN") {
          tobj = SubsetData(h1, ident.use=ccl.name)
        } else if(test.dir == "OUT") {
          tobj = SubsetData(h1, ident.remove=ccl.name)
        } else {
          stop("test.dir should be IN or OUT")
        }
      } # end if
      cat(ccl.name, " ")
      df.test = rbind(df.test,
                      TestGeneSig(tobj, sig.list, "sampleid", df.subj, 
                                  "response", c("low","high"), test.to.print = NA) %>% 
                        dplyr::mutate(cluster = ccl.name, test.dir, clustering)
      )
      cat("\n")
    } # ccl.name
  } # test.dir
} # clustering

fwrite(df.test, file.path(dn.out, "test_sig.genes_in_clusters_high_vs_low_responders.txt"), sep="\t")

