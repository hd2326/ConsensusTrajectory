DensityBasedManifoldDetection <- function(dist, radius, nobs, diag = F){
  dist <- as.matrix(dist)
  adjm <- matrix(0, nrow(dist), ncol(dist), dimnames = list(rownames(dist), colnames(dist)))
  for (i in 1:nrow(adjm)){
    if (sum(dist[i, ] < radius) > nobs) adjm[i, which(dist[i, ] < radius)] <- 1
  }
  if (diag) diag(adjm) <- 0
  return (adjm)
}
