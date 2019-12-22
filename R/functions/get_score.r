# Calculate score from genes x samples matrix as average z-score
get_score <- function(x) {
  x = t(scale(t(x)))
  return (colMeans(x, na.rm=T))
}
