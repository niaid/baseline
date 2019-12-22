library(data.table)

dn.out = file.path(PROJECT_DIR, "generated_data/brown-module-leading-genes-low-mid-high-eigengene")
dir.create(dn.out, showWarnings = F)

st = list({
     r = fread(file.path(PROJECT_DIR, "generated_data/brown-module-leading-genes-low-mid-high-eigengene/H1N1-zscore-AUC-low-mid-high-brown-leading-edge-wrt-CD38plusCellinFlu-87gene.csv"))
     r[,Study:="H1N1"]; r
}, {
     r = fread(file.path(PROJECT_DIR, "generated_data/brown-module-leading-genes-low-mid-high-eigengene/SDY212-zscore-AUC-low-mid-high-brown-leading-edge-wrt-CD38plusCellinFlu-87gene.csv"))
     r[,Study:="SDY212"]; r
}, {
     r = fread(file.path(PROJECT_DIR, "generated_data/brown-module-leading-genes-low-mid-high-eigengene/SDY400-zscore-AUC-low-mid-high-brown-leading-edge-wrt-CD38plusCellinFlu-87gene.csv"))
     r[,Study:="SDY400"]; r
}, {
     r = fread(file.path(PROJECT_DIR, "generated_data/brown-module-leading-genes-low-mid-high-eigengene/SDY404-zscore-AUC-low-mid-high-brown-leading-edge-wrt-CD38plusCellinFlu-87gene.csv"))
     r[,Study:="SDY404"]; r
})

st = rbindlist(st)
setnames(st,"MEbrown", "brown.L.E.eigengene")
setnames(st,"MEgrey", "grey.L.E.eigengene")
fn = file.path(dn.out, "four-flu-CD38plus-10gene-and-brown-mod-leading-edge-low-mid-high.csv")
fwrite(file=fn, st, quote=T)
st[,Study:=as.factor(Study)]

CD38sig.H1N1 = fread(file.path(PROJECT_DIR, "generated_data/CHI/CHI_cd38_ge_sig_score_day0.txt"))
setnames(CD38sig.H1N1,"subject","Sample")
setnames(CD38sig.H1N1,"score", "CD38_score")

CD38sig.HIPC = fread(file.path(PROJECT_DIR, "generated_data/HIPC/HIPC_cd38_ge_sig_score.txt"))[which(!(Response %in% "middle"))]
setnames(CD38sig.HIPC, "subject", "Sample")
stopifnot(!any(duplicated(CD38sig.HIPC$subject)))

tmp = rbindlist(list(CD38sig.H1N1[,c("Sample", "CD38_score"),wi=F], CD38sig.HIPC[,c("Sample", "CD38_score"),wi=F]))
st = merge(st, tmp, by="Sample")

saveRDS(st, file.path(dn.out, "st.rds"))

