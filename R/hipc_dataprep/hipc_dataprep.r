library(Biobase)

fn.esets = file.path(PROJECT_DIR, "generated_data","HIPC", "HIPC_IS_esets.rds")
eset_list = readRDS(fn.esets)

studies = names(eset_list)

for (sdy in studies) {
  eset = eset_list[[sdy]]
  dat = exprs(eset)
  info = pData(eset)
  
  si = info$Age.class == "young"
  dat = dat[,si]
  info = info[si,] %>% 
    dplyr::select(Study=study, sample=subject, Age, adjMFC_class = young_fc_res_max_d30) %>% 
    separate(sample, c("subject","time"), sep="_", remove = F) #%>% 

  if(sdy=="SDY404") { # updated data
    fn.t = file.path(PROJECT_DIR, "data", "HIPC", "phenotypes",
                     paste0(sdy, "_young_hai_titer_table_2016.txt"))
    titer = fread(fn.t) %>% dplyr::rename(subject=IDs)
    info = info %>% select(-adjMFC_class) %>% 
      left_join(titer %>% dplyr::select(subject, adjMFC_class=fc_res_max_d30), 
                by="subject") 
  }
  info = info %>% 
    mutate(Response = ifelse(adjMFC_class==0, "low",
                      ifelse(adjMFC_class==2, "high", "middle")))
  
  
  fn.ge = file.path(PROJECT_DIR, "generated_data", "HIPC",
                    sprintf("%s_GE_matrix.txt",sdy))
  dat %>% as.data.frame() %>% tibble::rownames_to_column("ID") %>%
    fwrite(fn.ge, sep="\t", quote=T)
  
  fn.si = file.path(PROJECT_DIR, "generated_data", "HIPC",
                    sprintf("%s_sample_info.txt",sdy))
  
  fwrite(info, fn.si, sep="\t", quote=T)
  
}
  