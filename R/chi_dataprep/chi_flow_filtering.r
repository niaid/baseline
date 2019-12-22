fn.poB = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_percent_of_B.txt")
flow.poB = fread(fn.poB) %>% tibble::column_to_rownames("sample") %>% 
  data.matrix()

fn.info = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_sample_info.txt")
flow.info = fread(fn.info) %>% 
  mutate(subject = as.character(subject))


i.flag = flow.info$Flag==0

fn.poB.f = sub(".txt", "_filtered.txt", fn.poB)
fn.info.f = sub(".txt", "_filtered.txt", fn.info)

flow.poB.f = flow.poB[i.flag,] %>% as.data.frame() %>% tibble::rownames_to_column("sample")
flow.info.f = flow.info[i.flag,]

fwrite(flow.poB.f, fn.poB.f, sep="\t", quote=T)
fwrite(flow.info.f, fn.info.f, sep="\t", quote=T)


i.d0 = flow.info.f$time == 0

fn.poB.f.d0 = sub(".txt", "_day0.txt", fn.poB.f)
fn.info.f.d0 = sub(".txt", "_day0.txt", fn.info.f)

flow.poB.f.d0 = flow.poB.f[i.d0,]
flow.info.f.d0 = flow.info.f[i.d0,]

fwrite(flow.poB.f.d0, fn.poB.f.d0, sep="\t", quote=T)
fwrite(flow.info.f.d0, fn.info.f.d0, sep="\t", quote=T)
