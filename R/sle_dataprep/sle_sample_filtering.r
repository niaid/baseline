fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene.txt")
dat = fread(fn.ge, data.table = F)

fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2.txt")
info = fread(fn.si)

si.sle = info$DISEASE == "SLE"

dat.sle = dat[,c(1,which(si.sle)+1)]
info.sle = info[si.sle,]

fn.ge.sle = sub(".txt", "_sle.txt", fn.ge)
fn.si.sle = sub(".txt", "_sle.txt", fn.si)

fwrite(dat.sle, fn.ge.sle, sep="\t", quote=T)
fwrite(info.sle, fn.si.sle, sep="\t", quote=T)

si.lowDA = info.sle$DA == "low"
dat.lowDA = dat.sle[,c(1,which(si.lowDA)+1)]
info.lowDA = info.sle[si.lowDA,]

fn.ge.lowDA = sub(".txt", "_lowDA.txt", fn.ge.sle)
fn.si.lowDA = sub(".txt", "_lowDA.txt", fn.si.sle)

fwrite(dat.lowDA, fn.ge.lowDA, sep="\t", quote=T)
fwrite(info.lowDA, fn.si.lowDA, sep="\t", quote=T)
