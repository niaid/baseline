# For genes with multiple probesets, pick probeset most correlated to the PC1 of all probesets.
# verified OK for random genes.
pick.probeset <- function(exprset,gene.map.file)
{
  library(Biobase)
  library(plyr)
  probinf = featureNames(exprset)
  dat = exprs(exprset)

  map = read.table(gene.map.file, sep="\t", header=T, stringsAsFactors=F)
  # overlap(map[,1], probinf)
  map = map[(map[,1] %in% probinf) & !is.na(map[,2]),]; 
        stopifnot(!duplicated(probinf), !duplicated(map[,1]), map[,1:2] != "", !is.na(map[,1:2]))
  map$idx = match(map[,1], probinf);  stopifnot(!is.na(map$idx))
  map = ddply(map, 2, function(df) 
        { 
          if (nrow(df)==1) {  
            return(data.frame(df, pc1varexpl=NA, nprobeset=1,
                              cor2pc1=NA, stringsAsFactors=F))
          }
          datgene = t(dat[df$idx,])
          pr = prcomp(datgene)
          if(all(pr$sdev==0)) {
            return(data.frame(df[1,], pc1varexpl=NA, nprobeset=nrow(df), 
                       cor2pc1=NA, stringsAsFactors=F))
          }
          pc1varexpl = 100*(pr$sdev[1])^2/sum(pr$sdev^2)
          cors = abs(cor(pr$x[,"PC1",drop=F], datgene))
          data.frame(df[which.max(cors),], pc1varexpl=pc1varexpl, nprobeset=nrow(df), cor2pc1=cors[which.max(cors)], stringsAsFactors=F)
        })
  stopifnot(!duplicated(map[,2]), !duplicated(map[,1]))
  write.table(map, file=sub('(\\.[^\\.]*)$','_PC1\\1',gene.map.file), sep="\t", quote=F, row.names=F)
  
  return(0)
} 
