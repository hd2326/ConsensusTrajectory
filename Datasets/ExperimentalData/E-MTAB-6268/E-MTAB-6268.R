library(igraph)
library(ape)
library(dbscan)
library(DDRTree)
source("../../../Algorithms/DensityBasedManifoldDetection/DensityBasedManifoldDetection.R")
source("../../ColorGradient.R")
load("E-MTAB-6268.rda")
den <- topo.colors(20)[sapply(colSums(as.matrix(dist(mds)) < 0.02), function(x) max(1, min(x, 20)))]
col <- rainbow(6)[as.factor(colnames(tpm))]

mst <- as.matrix(ape::mst(dist(mds)))
ddr <- as.matrix(DDRTree::DDRTree(t(mds), dimensions = 2)$stree) > 0
knn <- lapply(1:10, function(k, d){
  knn <- dbscan::kNN(dist(d), k = k, sort = F)$id
  knn <- cbind(as.numeric(sapply(1:nrow(knn), function(x) rep(x, k))), as.numeric(t(knn)))
  as.matrix(igraph::as_adj(graph_from_edgelist(knn)))
}, d=mds)
dbmd <- DensityBasedManifoldDetection(dist(mds), 0.1, 10)
dbmd <- dbmd + t(dbmd) == 2
#algorithms

consensus <- mst + ddr + Reduce("+", knn) + dbmd
svd <- svd(consensus)
D <- diag(svd$d)
diag(D)[30:ncol(D)] <- 0
consensus1 <- svd$u %*% D %*% t(svd$v)
consensus1 <- consensus1 >= 1
consensus1 <- igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected")
modularity <- igraph::cluster_walktrap(consensus1, weights = NULL)
centrality <- igraph::centralization.betweenness(consensus1)
startpoints <- sapply(1:max(modularity$membership), function(i, centrality, modularity, g){
  c <- structure(centrality$res, names = 1:igraph::vcount(g))
  names(c[modularity$membership == i])[which.min(c[modularity$membership == i])]
}, centrality=centrality, modularity=modularity, g=consensus1)
distance <- lapply(startpoints, function(p, g){
  sapply(1:igraph::vcount(g), function(i, s, g) igraph::distances(g, s, i), g=g, s=as.numeric(p))
}, g=consensus1)

pdf("Fig2.pdf", width = 10, height = 5)
par(mfrow = c(2, 4), mar = c(2, 2, 2, 2))
set.seed(1)
plot(consensus1, vertex.color = col, vertex.label = "", vertex.size = 6, main = "E-MTAB-6268")
legend("topleft", legend = levels(as.factor(colnames(tpm)))[c(1, 3, 5, 2, 4)], col = rainbow(6)[c(1, 3, 5, 2, 4)], pch = 16, pt.cex = 2, bty = "n", border = NA)

set.seed(1)
plot(consensus1, vertex.color = topo.colors(max(modularity$membership))[modularity$membership], vertex.label = "", vertex.size = 6, main = "Clustering Analysis")
legend("topleft", legend = paste("C", 1:max(modularity$membership)), col = topo.colors(max(modularity$membership)), pch = 16, pt.cex = 2, bty = "n", border = NA)

size <- log2(centrality$res+1)
size[c(as.numeric(startpoints[[1]]), as.numeric(startpoints[[2]]), as.numeric(startpoints[[5]]))] <- max(size)+1
shape <- rep("circle", igraph::vcount(consensus1))
shape[as.numeric(startpoints[[5]])] <- "square"
col <- heat.colors(max(distance[[5]])+1, alpha = 0.8)[distance[[5]]+1]
col[c(as.numeric(startpoints[[1]]), as.numeric(startpoints[[2]]), as.numeric(startpoints[[5]]))] <- "Blue"
set.seed(1)
plot(consensus1, vertex.color = col, vertex.label = "", vertex.size = size, vertex.shape = shape, main = "Lineage Analysis")
legend("bottomleft", legend = c("High Betweenness", "Low Betweenness"), pch = 16, pt.cex = c(2, 1), bty = "n", border = NA)
legend("topleft", legend = c("Start Cell", "End Cells", "Early Pseudotime", "Late Pseudotime"), col = c("Blue", "Blue", "Red", "Yellow"), pt.cex = 2, pch = c(15, 16, 16, 16), bty = "n", border = NA)

col <- rep(1, igraph::vcount(consensus1))
col[modularity$membership == 4] <- "Blue"
size <- rep(6, igraph::vcount(consensus1))
size[modularity$membership == 4] <- 15
shape <- rep("circle", igraph::vcount(consensus1))
shape[modularity$membership == 4] <- "square"
set.seed(1)
plot(consensus1, vertex.color = col, vertex.label = "", vertex.size = size, vertex.shape = shape, main = "Branch Analysis")
legend("bottomleft", legend = c("Branching Cells"), pch = 15, pt.cex = 2, col = "Blue", bty = "n", border = NA)

#sort(rowMeans(tpm[, branch == 0]) - rowMeans(tpm[, branch > 0]), decreasing = T)[1:10]
set.seed(1)
plot(consensus1, vertex.color = expColor("HAPLN1", tpm), vertex.label = "", vertex.size = 6, main = "HAPLN1 (#1 for Branch Cells)")
legend("topleft", legend = c("High Expression", "Low Expression"), pch = 16, col = c("Red", "Grey"), pt.cex = 2, bty = "n", border = NA)
set.seed(1)
plot(consensus1, vertex.color = expColor("TMEM88", 1.3^tpm), vertex.label = "", vertex.size = 6, main = "TMEM88 (#2 for Branch Cells)")
legend("topleft", legend = c("High Expression", "Low Expression"), pch = 16, col = c("Red", "Grey"), pt.cex = 2, bty = "n", border = NA)
set.seed(1)
plot(consensus1, vertex.color = expColor("COL1A1", 1.3^tpm), vertex.label = "", vertex.size = 6, main = "COL1A1 (Fibroblast)")
legend("topleft", legend = c("High Expression", "Low Expression"), pch = 16, col = c("Red", "Grey"), pt.cex = 2, bty = "n", border = NA)
set.seed(1)
plot(consensus1, vertex.color = expColor("TNNT2", 1.3^tpm), vertex.label = "", vertex.size = 6, main = "TNNT2 (Cardiomyocyte)")
legend("topleft", legend = c("High Expression", "Low Expression"), pch = 16, col = c("Red", "Grey"), pt.cex = 2, bty = "n", border = NA)
dev.off()
#Fig2

pdf("SupFig5.pdf", width = 10, height = 5)
par(mfcol = c(2, 4))
consensus <- mst + ddr + Reduce("+", knn) + dbmd
svd <- svd(consensus)
plot(mds, col = col, xlab = "Dim1", ylab = "Dim2", main = "E-MTAB-6268")
legend("topleft", legend = levels(as.factor(colnames(tpm))), fill = rainbow(6), bty = "n", border = NA)
plot(svd$d, xlab = "Rank", ylab = "SV", xaxt = "n", main = "SVD")
axis(side = 1, at = c(5, 30, 100), las = 2)
par(mar = c(1, 1, 3, 1))
for (i in c(5, 30, 100)){
  D <- diag(svd$d)
  diag(D)[i:ncol(D)] <- 0
  consensus1 <- svd$u %*% D %*% t(svd$v)
  consensus1 <- consensus1 >= 1
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected"), vertex.color = col, vertex.label = "", vertex.size = 6, main = paste("First", i, "SVs"))
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected"), vertex.color = den, vertex.label = "", vertex.size = 6, main = "Local Density")
}
dev.off()
#SupFig5
