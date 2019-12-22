library(Biobase)

fn.esets = file.path(PROJECT_DIR, "generated_data","HIPC", "HIPC_IS_esets.rds")
eset_list = readRDS(fn.esets)

studies = names(eset_list)

for (sdy in studies) {
  
  eset = eset_list[[sdy]]
  
  # generate probe to gene map
  probes.map = featureData(eset)@data %>% 
    tibble::rownames_to_column("ID") %>% 
    dplyr::filter(gene_symbol != "")
  colnames(probes.map) = c("ID","gene")
  fn.map = file.path(PROJECT_DIR, "generated_data", "HIPC",
                     paste0(sdy,"_probe_map.txt"))
  fwrite(probes.map, file=fn.map, sep="\t", quote=T)
  
  source(file.path(PROJECT_DIR, "R/functions/pick.probeset.r"))
  pick.probeset(eset,fn.map)
}
