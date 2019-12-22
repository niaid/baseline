library(data.table)

brown.LE.isv050 = fread(file.path(PROJECT_DIR, "generated_data/fgsea_with_wgcna_modules/brown-leading-edge-isv050-87genes.txt"))$gene

#######
###
sdy212         = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY212_GE_matrix_gene_day0.txt"))
sdy212.samples = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY212_sample_info_day0.txt"))
setnames(sdy212.samples, "subject", "SubjectID")
sdy212.d0 = sdy212[,c(1, grep("_d0$", colnames(sdy212))), wi=F]
setnames(sdy212.d0, 1:ncol(sdy212.d0), gsub("_d0", "", colnames(sdy212.d0)))

sdy400         = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY400_GE_matrix_gene_day0.txt"))
sdy400.samples = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY400_sample_info.txt"))
setnames(sdy400.samples, "subject", "SubjectID")
sdy400.d0 = sdy400[,c(1, grep("_d0$", colnames(sdy400))), wi=F]
setnames(sdy400.d0, 1:ncol(sdy400.d0), gsub("_d0", "", colnames(sdy400.d0)))

sdy404         = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY404_GE_matrix_gene_day0.txt"))
sdy404.samples = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY404_sample_info_day0.txt"))
setnames(sdy404.samples, "subject", "SubjectID")
sdy404.d0 = sdy404[,c(1, grep("_d0$", colnames(sdy404))), wi=F]
setnames(sdy404.d0, 1:ncol(sdy404.d0), gsub("_d0", "", colnames(sdy404.d0)))

H1N1         = fread(file.path(PROJECT_DIR, "generated_data/CHI/CHI_GE_matrix_gene.txt"))
H1N1.samples = fread(file.path(PROJECT_DIR, "generated_data/CHI/CHI_sample_info_2_CD38hi.txt"))
setnames(H1N1.samples, "subject", "subject.id")
H1N1.samples[,SubjectID:=as.character(subject.id)]
H1N1.d0 = H1N1[,c(1, grep("_day0$", colnames(H1N1))), wi=F]
setnames(H1N1.d0, 1:ncol(H1N1.d0), gsub("_day0", "", colnames(H1N1.d0)))


