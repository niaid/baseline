library(Biobase)

fn.in = "data/CHI/expression/corrected.batch.txt"
dat = fread(file.path(PROJECT_DIR, fn.in), header=T, data.table=F)

info = dat[,c(1,3)]

info = info %>% 
  mutate(subject = as.character(subject)) %>% 
  mutate(time.point = paste0(ifelse(time<0,"pre","day"),abs(time))) %>% 
  mutate(sample = paste(subject, time.point, sep="_"))

dat2 = dat[,-c(1:3)] %>% data.matrix() 
rownames(dat2) = info$sample
dat2 = t(dat2)

# output gene expression data and sample info
dir.create("generated_data/CHI", showWarning = F)
fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix.txt")
fn.si = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_sample_info.txt")
dat2 %>% as.data.frame() %>% tibble::rownames_to_column("ID") %>% 
  fwrite(fn.ge, sep="\t", quote=T)
fwrite(info, fn.si, sep="\t", quote=T)

