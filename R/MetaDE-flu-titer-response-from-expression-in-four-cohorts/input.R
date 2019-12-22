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
dat1 = fread(file.path(PROJECT_DIR, "generated_data/CHI/CHI_GE_matrix_gene.txt"))
dat2 = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY212_GE_matrix_gene_day0.txt"))
dat3 = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY400_GE_matrix_gene_day0.txt"))
dat4 = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY404_GE_matrix_gene_day0.txt"))

## convert to upper case
dat1[,gene:=toupper(gene)]; setkey(dat1, gene); dim(dat1)
dat2[,gene:=toupper(gene)]; setkey(dat2, gene); dim(dat2)
dat3[,gene:=toupper(gene)]; setkey(dat3, gene); dim(dat3)
dat4[,gene:=toupper(gene)]; setkey(dat4, gene); dim(dat4)

## find duplicated symbols and 
dup1 = which(duplicated(dat1[[1]]))
dup2 = which(duplicated(dat2[[1]]))
dup3 = which(duplicated(dat3[[1]]))
dup4 = which(duplicated(dat4[[1]]))

# dat2[[1]][dup2]# %in% WGCNA[[1]]
# dat3[[1]][dup3]# %in% WGCNA[[1]]
# dat4[[1]][dup4]# %in% WGCNA[[1]]

# only two duplicated symbols
# "RG9MTD1"  # duplicated symbol within dat2, dat3, and dat4, but not WGCNA
# "GCOM1"    # duplicated symbol within       dat3, and dat4, and in grey module from WGCNA

clean.symbol = function(dat,g) {
  #g = "GCOM1"
  #dat = dat3
  idx = which(dat[[1]] %in% g)
  if(length(idx)>1) {
    idx.to.keep = which.max( diag(var(t(dat[idx,-1,wi=F]))))
    dat = dat[-idx[-idx.to.keep]]
  }
  return(dat)
}
dim(dat3)
dat3 =clean.symbol(dat3, "GCOM1")
dat4 =clean.symbol(dat4, "GCOM1")
dim(dat3)
dim(dat4)

## for duplicated symbols, keep the one with bigger variance
common.genes = Reduce(intersect, list(dat1$gene, dat2$gene, dat3$gene, dat4$gene), dat1$gene)
common.genes = setdiff(common.genes, "RG9MTD1")
summary(common.genes)

fn = file.path(dn.fig, "Number-of-genes-in-4flu-datasets.pdf")
pdf(fn)
       par(mar=c(7,15,3,5))
       tmp = c(
               length(intersect(common.genes, WGCNA[[1]])),
               nrow(WGCNA),
               length(common.genes), 
               nrow(dat1),
               nrow(dat2),
               nrow(dat3),
               nrow(dat4))
       names(tmp) = c("Overlap between four flu datasets\nand SLE low-DA datasets",
                      "SLE low-DA datasets",
                      "Common to the four\ndatasets above",
                      "H1N1",
                      "SDY212",
                      "SDY400",
                      "SDY404")
       main="Number of genes in each dataset"
       barplot(tmp, horiz=T, las=2, main=main, xlab="# genes")
       text(tmp,((1:7)-0.5)*1.2,tmp,pos=2)
dev.off()


if(1) {
  cnt = sapply(unique(WGCNA$Module), function(x) { 
           y=WGCNA[which(Module %in% x)]$Symbol; 
           c(length(y),(length(intersect(common.genes, y))))})
  cnt = cnt[, colnames(cnt) != "grey"]
  
  pdf(file.path(dn.fig, "WGCNA-genes.pdf"))
         par(mar=c(5,8,5,5))
         o = order(-cnt[1,])
         n = ncol(cnt)
         main = "Number of genes in WGCNA modules"
         barplot(cnt[1,o],horiz=T,col=colnames(cnt)[o], log="x", las=2, xlim=c(20,1000), main=main)
         text(cnt[1,o], ((1:n)-0.5)*1.2, cnt[1,o], pos=4)
         text(cnt[1,o[1]]-1000,0.6,cnt[1,o[1]], pos=2)
  dev.off()

  pdf(file.path(dn.fig, "WGCNA-genes-in-four-flu-datasets.pdf"))
         par(mar=c(5,8,5,5))
         o = order(-cnt[1,])
         n = ncol(cnt)
         main = "Genes common to flu datasets also found in WGCNA modules"
         barplot(cnt[2,o],horiz=T,col=colnames(cnt)[o], log="x", las=2, xlim=c(3,1000), main=main)
         text(cnt[2,o][-1], (((1:n)-0.5)*1.2)[-1], cnt[2,o][-1], pos=4)
         text(cnt[2,o[1]]-1000,0.6,cnt[2,o[1]], pos=2)
  dev.off()
}

if(1) {
  annot1 = fread(file.path(PROJECT_DIR, "generated_data/CHI/CHI_sample_info_2_CD38hi.txt"))[which(Response %in% c("low","high"))][which(time %in% 0)]
  label1 = annot1[,c("sample","Response"),wi=F]
  label1[which(Response %in% "low"), class:=0]
  label1[which(Response %in% "high"), class:=1]
  setkey(label1,class)  # use label1$class in MetaDE
  
  exp1 = as.matrix(dat1[common.genes][,label1$sample,wi=F])
  rownames(exp1) <- common.genes
  dim(exp1)
}

if(1) {
  annot2 = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY212_sample_info_day0.txt"))
  setnames(annot2, "subject", "SubjectID")
  label2 = annot2[which(Response %in% c("low","high")), c("SubjectID","Response"),wi=F]
  label2[which(Response %in% "low"), class:=0]
  label2[which(Response %in% "high"), class:=1]
  setkey(label2, class)

  exp2 = as.matrix(dat2[common.genes][,paste0(label2$SubjectID,"_d0"),wi=F])
  rownames(exp2) <- common.genes
}

if(1) {
  annot3 = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY400_sample_info_day0.txt"))
  setnames(annot3, "subject", "SubjectID")
  label3 = annot3[which(time %in% "d0" & Response %in% c("low","high")), c("SubjectID","Response"),wi=F]
  label3[which(Response %in% "low"), class:=0]
  label3[which(Response %in% "high"), class:=1]
  setkey(label3, class)

  exp3 = as.matrix(dat3[common.genes][,paste0(label3$SubjectID,"_d0"),wi=F])
  rownames(exp3) <- common.genes
}


if(1) {
  annot4 = fread(file.path(PROJECT_DIR, "generated_data/HIPC/SDY404_sample_info_day0.txt"))
  setnames(annot4, "subject", "SubjectID")
  label4 = annot4[which(time %in% "d0" & Response %in% c("low","high")), c("SubjectID","Response"),wi=F]
  label4[which(Response %in% "low"), class:=0]
  label4[which(Response %in% "high"), class:=1]
  setkey(label4, class)

  exp4 = as.matrix(dat4[common.genes][,paste0(label4$SubjectID,"_d0"),wi=F])
  rownames(exp4) <- common.genes
}

if(1) {
  pdf(file.path(dn.fig, "Number-of-subjects-in-4flu-datasets.pdf"))
    tmp = c(table(label1$class),
           table(label2$class),
           table(label3$class),
           table(label4$class)
           )
    names(tmp) = rep(c("Low", "High"), 4)
    names(tmp) = paste(c(rep("H1N1",2),
                        rep("SDY212", 2),
                        rep("SDY400", 2),
                        rep("SDY404", 2)), names(tmp))
    par(mar=c(5,10,5,2))
    # plot.new()
    main = sprintf("Number of subjects in the four flu datasets")
    barplot(tmp, horiz=T, las=1, xlim=c(0,16), main=main, xlab="#Subject")
    text(tmp, ((1:8)-0.5)*1.2, tmp, pos=4)
  dev.off()
}


