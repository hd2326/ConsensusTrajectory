library(igraph)
library(ape)
library(dbscan)
library(DDRTree)
source("../../../Algorithms/DensityBasedManifoldDetection/DensityBasedManifoldDetection.R")
load("E-MTAB-3929.rda")
set.seed(1)
ind <- sapply(unique(anno), function(x, anno) sample(which(anno == x), 50, replace = F), anno=anno)
col <- rainbow(6)[as.factor(anno[ind])]
den <- topo.colors(20)[sapply(colSums(as.matrix(dist(pca$x[, 1:2])) < 1), function(x) max(1, min(x, 20)))]

mst <- as.matrix(ape::mst(dist(pca$x[ind, 1:2])))
ddr <- as.matrix(DDRTree::DDRTree(t(pca$x[ind, 1:2]), dimensions = 2)$stree) > 0
knn <- lapply(1:10, function(k, d){
  knn <- dbscan::kNN(dist(d), k = k, sort = F)$id
  knn <- cbind(as.numeric(sapply(1:nrow(knn), function(x) rep(x, k))), as.numeric(t(knn)))
  as.matrix(igraph::as_adj(graph_from_edgelist(knn)))
}, d=pca$x[ind, 1:2])
dbmd <- DensityBasedManifoldDetection(dist(pca$x[ind, 1:2]), 10, 10)
dbmd <- dbmd + t(dbmd) == 2
#algorithms

pdf("SupFig6.pdf", width = 12, height = 6)
par(mfcol = c(2, 4))
consensus <- mst + ddr + Reduce("+", knn) + dbmd
svd <- svd(consensus)
plot(pca$x[ind, 1:2], col = col, xlab = "Dim1", ylab = "Dim2", main = "E-MTAB-3929")
legend("topleft", legend = levels(as.factor(anno)), fill = rainbow(6), bty = "n", border = NA)
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
  plot(igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected"), vertex.color = den, vertex.label = "", vertex.size = 6, main = "Local Density")
}
dev.off()
#SupFig6

D <- diag(svd$d)
diag(D)[30:ncol(D)] <- 0
consensus1 <- svd$u %*% D %*% t(svd$v)
consensus1 <- consensus1 >= 1
consensus1 <- igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected")
modularity <- igraph::cluster_walktrap(consensus1, weights = NULL)
centrality <- igraph::centralization.betweenness(consensus1)
distance <- sapply(1:igraph::vcount(consensus1), function(i, s, g) igraph::distances(g, s, i), g=consensus1, s=which.min(centrality$res))
pdf("SupFig7.pdf", width = 9, height = 6)
par(mfrow = c(2, 3))
set.seed(1)
size <- log2(centrality$res)
size[which.min(size)] <- max(size)+1
shape <- c("circle", "square")[(centrality$res == min(centrality$res))+1]
plot(consensus1, vertex.color = heat.colors(max(distance)+2)[distance+1], vertex.label = "", vertex.size = size, vertex.shape = shape, main = "Lineage Analysis")
legend("bottomright", legend = c("High Betweenness", "Low Betweenness"), pch = 16, pt.cex = c(2, 1), bty = "n", border = NA)
legend("topleft", legend = c("Start Cell", "Early Pseudotime","Late Pseudotime"), pch = c(15, 16, 16), col = c("Red", "Red", "Yellow"), pt.cex = 2, bty = "n", border = NA)
plot(distance, rpkm["GATA3", ind], xlab = "Distance", ylab = "Expression", col = rainbow(6)[as.factor(anno[ind])], pch = 16, main = "GATA3")
plot(distance, rpkm["DPPA5", ind], xlab = "Distance", ylab = "Expression", col = rainbow(6)[as.factor(anno[ind])], pch = 16, main = "DPPA5")
legend("bottomright", legend = levels(as.factor(anno)), col = rainbow(6), pch = 16, pt.cex = 2, bty = "n", border = NA)

set.seed(1)
plot(consensus1, vertex.color = topo.colors(max(modularity$membership))[modularity$membership], vertex.label = "", vertex.size = 6, main = "Clustering Analysis")
legend("bottomright", legend = paste("C", 1:max(modularity$membership)), col = topo.colors(max(modularity$membership)), pch = 16, pt.cex = 2, bty = "n", border = NA)
set.seed(1)
plot(consensus1, vertex.color = col, vertex.label = "", vertex.size = 6, main = "")
legend("topleft", legend = levels(as.factor(anno)), col = rainbow(6), pch = 16, pt.cex = 2, bty = "n", border = NA)
boxplot(list(rpkm["GATA3", ind][modularity$membership == 1], rpkm["GATA3", ind][modularity$membership == 2],
             rpkm["GATA3", ind][modularity$membership == 3], rpkm["GATA3", ind][modularity$membership == 4],
             rpkm["DPPA5", ind][modularity$membership == 1], rpkm["DPPA5", ind][modularity$membership == 2],
             rpkm["DPPA5", ind][modularity$membership == 3], rpkm["DPPA5", ind][modularity$membership == 4]),
        border = rep(topo.colors(max(modularity$membership)), 2), names = NA, ylab = "Expression",xaxt = "n", main = "")
axis(side = 1, at = c(2.5, 6.5), labels = c("GATA3", "DPPR5"), tick = F)
dev.off()
#SupFig7
