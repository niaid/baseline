fn.ge = file.path(PROJECT_DIR, "generated_data", "YF", "YF_GE_matrix.txt")
dat = fread(fn.ge, sep="\t")

fn.map.pc1 = file.path(PROJECT_DIR, "generated_data", "YF", "YF_probe_map_PC1.txt")
probes.map.pc1 = fread(fn.map.pc1, sep="\t")
gi = match(probes.map.pc1$ID, dat$ID)
dat.gene = dat[gi,]
dat.gene$ID  <- probes.map.pc1$gene
names(dat.gene)[1] = "gene"

fn.ge.gene = sub(".txt", "_gene.txt", fn.ge)
fwrite(dat.gene, fn.ge.gene, sep="\t", quote=T)
