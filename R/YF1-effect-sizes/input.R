library(data.table)

dn.out = file.path(PROJECT_DIR, "generated_data/MetaDE")
dir.create(dn.out, showWarnings = F)
dn.fig = file.path(PROJECT_DIR, "figure_generation/MetaDE")
dir.create(dn.fig, showWarnings = F)

WGCNA = fread(file.path(PROJECT_DIR, "generated_data/WGCNA-modules-from-SLE-low-DA/SLE-low-34sbj-9601probes-gene.in.module.minModSize20.signed-hybrid.txt"))
WGCNA[,Symbol:=toupper(Symbol)]

### STEP 1. read data ###
###         MetaDE assumptions:
###         1. each column of expression matrix is a sample 
###         2. labels are numeric vectors
dat0 = fread(file.path(PROJECT_DIR, "generated_data/YF/YF_GE_matrix_gene_day0.txt"))
dat0[,gene:=toupper(gene)]
setkey(dat0, gene)
setnames(dat0, colnames(dat0), gsub(".CEL$", "", colnames(dat0)))
stopifnot(!any(duplicated(colnames(dat0))))

annot0 = fread(file.path(PROJECT_DIR, "generated_data/YF/YF_sample_info_day0.txt"))
annot0[,Trial:=as.character(Trial)]
annot1 = annot0[which(Trial=="1" & Time==0 & Response %in% c("low","high"))]
annot2 = annot0[which(Trial=="2" & Time==0 & Response %in% c("low","high"))]
annot1[which(Response %in% "low"), class:=0]
annot1[which(Response %in% "high"), class:=1]
annot2[which(Response %in% "low"), class:=0]
annot2[which(Response %in% "high"), class:=1]

pdf(file.path(dn.fig, "YF_sample-overview.pdf"))
       tmp = c(table(annot1$class),
               table(annot2$class))
       names(tmp) = c("Trial1 Low", "Trial1 High",
                      "Trial2 Low", "Trial2 High")
       par(mar=c(5,10,5,5))
       main="Yellow Fever Vaccine Trial Day 0 Samples\n(with high or low response)"
       barplot(tmp, horiz=T, las=1, main=main)
       text(tmp, ((1:4)-0.5)*1.2, tmp, pos=2)
dev.off()

pdf(file.path(dn.fig, "YF_sample-overview-all-timepoints.pdf"))
       tmp = annot0[,list(Sample=length(subject)),by=c("Trial","Time")]
       tmp[,Time:=factor(Time)]
       p = ggplot(tmp, aes(x=Trial, y=Sample, fill=Time)) + geom_bar(stat="identity", position="dodge") +
         scale_fill_discrete(name="Time") + labs(title="Yellow Fever Vaccine Dataset")
       print(p)
dev.off()

setkey(annot1, class, subject)
setkey(annot2, class, subject)

## extract trial 1 and 2 separately
dat1 = dat0[,c("gene", annot1$geo_accession),wi=F]
dat2 = dat0[,c("gene", annot2$geo_accession),wi=F]

## find duplicated symbols and 
dup1 = which(duplicated(dat1[[1]]))
dup2 = which(duplicated(dat2[[1]]))
# No duplicated symbols

common.genes = dat1$gene

if(1) {
  cnt = sapply(unique(WGCNA$Module), function(x) { 
           y=WGCNA[which(Module %in% x)]$Symbol; 
           c(length(y),(length(intersect(common.genes, y))))})
  cnt = cnt[, colnames(cnt) != "grey"]
  
  pdf(file.path(dn.fig, "WGCNA-genes-in-yellow-fever-datasets.pdf"))
         par(mar=c(5,8,5,5))
         o = order(-cnt[1,])
         n = ncol(cnt)
         main = "Genes common to yellow fever datasets also\nfound in WGCNA modules"
         barplot(cnt[2,o],horiz=T,col=colnames(cnt)[o], log="x", las=2, xlim=c(3,1000), main=main)
         text(cnt[2,o][-1], (((1:n)-0.5)*1.2)[-1], cnt[2,o][-1], pos=4)
         text(cnt[2,o[1]]-1000,0.6,cnt[2,o[1]], pos=2)
  dev.off()
}

if(1) {
  label1 = annot1[T] # use label1$class in MetaDE
  
  setkey(dat1,gene)
  exp1 = as.matrix(dat1[common.genes,-1,wi=F])
  rownames(exp1) <- common.genes
}

if(1) {
  label2 = annot2[T]

  setkey(dat2,gene)
  exp2 = as.matrix(dat2[common.genes,-1,wi=F])
  rownames(exp2) <- common.genes
}
