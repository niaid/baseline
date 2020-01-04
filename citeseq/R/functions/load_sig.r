load_sig <- function(data=NULL, col=NULL, desc.order=T, ntop) {
  if(is.null(data)) stop("Please enter file name")
  if(is.null(col)) stop("Please enter a column name to rank genes")
  if(is.character(data)) {
    if(file.exists(data)) {
      df.sig = read.table(data, sep="\t", header=T, row.names=1, stringsAsFactors=F)
    } else {
      stop(sprintf("File not found - %s", data))
    }
  } else if (is.data.frame(data)){
    df.sig = data
  } else {
    stop("Unsuported data. Specify a file name or a data frame")
  }
  if (is.character(col) & length(col)==1) {
    predictor = df.sig[,names(df.sig) == col]
  } else if (is.numeric(col) & length(col)==1) {
    predictor = df.sig[,col]
  } else {
    stop("Unsupported column variable. Use single column name or number")
  }
  names(predictor) = rownames(df.sig)
  predictor = sort(predictor, decreasing = desc.order)
  return(names(predictor)[1:ntop])
}
