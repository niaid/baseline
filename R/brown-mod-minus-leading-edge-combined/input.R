library(data.table)

H1N1 = fread(file.path(PROJECT_DIR, "generated_data/brown-mod-minus-leading-edge/H1N1-AUC-brown-minus-leading-edge-wrt-CD38plusCellinFlu-283gene.csv"))
SDY212 = fread(file.path(PROJECT_DIR, "generated_data/brown-mod-minus-leading-edge/SDY212-AUC-brown-minus-leading-edge-wrt-CD38plusCellinFlu-283gene.csv"))
SDY404 = fread(file.path(PROJECT_DIR, "generated_data/brown-mod-minus-leading-edge/SDY404-AUC-brown-minus-leading-edge-wrt-CD38plusCellinFlu-283gene.csv"))
SDY400 = fread(file.path(PROJECT_DIR, "generated_data/brown-mod-minus-leading-edge/SDY400-AUC-brown-minus-leading-edge-wrt-CD38plusCellinFlu-283gene.csv"))

H1N1[,Study:="NIH_2008"]
SDY212[,Study:="Stanford_2008"]
SDY404[,Study:="Yale_2011"]
SDY400[,Study:="Yale_2012"]

flu = rbindlist(list(H1N1, SDY212, SDY404, SDY400))
flu[,MEgrey:=NULL]
setnames(flu, "MEbrown", "brown.minus.LE")

CD38sig.H1N1 = fread(file.path(PROJECT_DIR, "generated_data/CHI/CHI_cd38_ge_sig_score_day0.txt"))
setnames(CD38sig.H1N1,"subject","Sample")
setnames(CD38sig.H1N1,"score", "CD38_score")

CD38sig.HIPC = fread(file.path(PROJECT_DIR, "generated_data/HIPC/HIPC_cd38_ge_sig_score.txt"))[which(!(Response %in% "middle"))]
setnames(CD38sig.HIPC, "subject", "Sample")
stopifnot(!any(duplicated(CD38sig.HIPC$subject)))

tmp = rbindlist(list(CD38sig.H1N1[,c("Sample", "CD38_score"),wi=F], CD38sig.HIPC[,c("Sample", "CD38_score"),wi=F]))
flu = merge(flu, tmp, by="Sample")

