library(igraph)
library(ape)
library(dbscan)
library(DDRTree)
source("../../../Algorithms/DensityBasedManifoldDetection/DensityBasedManifoldDetection.R")
load("E-MTAB-2805.rda")
set.seed(1)
ind <- sapply(unique(anno), function(x, anno) sample(which(anno == x), 50, replace = F), anno=anno)
col <- rainbow(4)[as.factor(anno[ind])]
den <- topo.colors(20)[sapply(colSums(as.matrix(dist(pca$x[, 1:2])) < 10), function(x) max(1, min(x, 20)))]

mst <- as.matrix(ape::mst(dist(pca$x[ind, 1:2])))
ddr <- as.matrix(DDRTree::DDRTree(t(pca$x[ind, 1:2]), dimensions = 2)$stree) > 0
knn <- lapply(1:10, function(k, d){
  knn <- dbscan::kNN(dist(d), k = k, sort = F)$id
  knn <- cbind(as.numeric(sapply(1:nrow(knn), function(x) rep(x, k))), as.numeric(t(knn)))
  as.matrix(igraph::as_adj(graph_from_edgelist(knn)))
}, d=pca$x[ind, 1:2])
dbmd <- DensityBasedManifoldDetection(dist(pca$x[ind, 1:2]), 20, 10)
dbmd <- dbmd + t(dbmd) == 2
#algorithms

pdf("SupFig8.pdf", width = 12, height = 6)
par(mfcol = c(2, 4))
consensus <- mst + ddr + Reduce("+", knn) + dbmd
svd <- svd(consensus)
plot(pca$x[ind, 1:2], col = col, xlab = "Dim1", ylab = "Dim2", main = "E-MTAB-2805")
legend("topleft", legend = levels(as.factor(anno)), fill = rainbow(4), bty = "n", border = NA)
plot(svd$d, xlab = "Rank", ylab = "SV", xaxt = "n", main = "SVD")
axis(side = 1, at = c(5, 30, 100))
par(mar = c(1, 1, 3, 1))
for (i in c(5, 30, 100)){
  D <- diag(svd$d)
  diag(D)[i:ncol(D)] <- 0
  consensus1 <- svd$u %*% D %*% t(svd$v)
  consensus1 <- consensus1 >= 1
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected"), vertex.color = col, vertex.label = "", vertex.size = 6, main = paste("First", i, "SVs"))
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected"), vertex.color = den, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "Local Density")
}
dev.off()
#SupFig8
