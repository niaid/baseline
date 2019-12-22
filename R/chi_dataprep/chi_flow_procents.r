fn.flow = file.path(PROJECT_DIR, "data", "CHI", "flow", "bp_flow_data.txt")
flow = fread(fn.flow) %>% tibble::column_to_rownames("Sample") %>% 
  data.matrix()

flow.poB = flow[,-(1:5)] / flow[,"B_cells"] * 100

flow.info = sub("^\\d+: (.+)\\.fcs","\\1",rownames(flow.poB)) %>% 
            gsub("pre ","-",.) %>% data.frame(sample=.) %>% 
            separate(sample, c("subject","time"), sep="_", remove = F) %>% 
            mutate(time = as.integer(time)) %>% 
            mutate(time.point = paste0(ifelse(time<0,"pre","day"),abs(time))) %>% 
            mutate(sample = paste(subject, time.point, sep="_"))

rownames(flow.poB) = flow.info$sample

fn.flag = file.path(PROJECT_DIR, "data", "CHI", "flow", "bp_flow_sample_flagged.txt")
flow.flag = fread(fn.flag)

flow.info = flow.info %>% left_join(flow.flag, by="sample")

fn.flow.poB = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_percent_of_B.txt")
fn.flow.si = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_sample_info.txt")

flow.poB %>% as.data.frame() %>% tibble::rownames_to_column("sample") %>% 
  fwrite(fn.flow.poB, sep="\t", quote=T)

fwrite(flow.info, fn.flow.si, sep="\t", quote=T)
