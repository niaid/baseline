fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix.txt")

dat = fread(fn.ge)

fn.map.pc1 = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_probe_map_PC1.txt")
probes.map.pc1 = read_tsv(fn.map.pc1)
gi = match(probes.map.pc1$ID,dat$ID)
dat.gene = dat[gi,]
dat.gene$ID  <- probes.map.pc1$gene
names(dat.gene)[1] = "gene"

fn.ge.gene = sub(".txt", "_gene.txt", fn.ge)
fwrite(dat.gene, fn.ge.gene, sep="\t", quote=T)
