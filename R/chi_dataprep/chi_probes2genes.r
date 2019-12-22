fn.ge = file.path(PROJECT_DIR, "generated_data", "CHI", "CHI_GE_matrix.txt")

dat = read_tsv(fn.ge)

fn.map.pc1 = file.path(PROJECT_DIR, "data", "CHI", "expression", "affy_hugene_1.0_ID_unique_PC1.txt")
probes.map.pc1 = read_tsv(fn.map.pc1)
gi = match(probes.map.pc1$probeid,dat$ID)
dat.gene = dat[gi,]
dat.gene$ID  <- probes.map.pc1$symbol
names(dat.gene)[1] = "gene"

fn.ge.gene = sub(".txt", "_gene.txt", fn.ge)
fwrite(dat.gene, fn.ge.gene, sep="\t", quote=T)
