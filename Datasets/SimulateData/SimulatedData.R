library(igraph)
library(ape)
library(dbscan)
library(DDRTree)
library(destiny)
source("../../Algorithms/DensityBasedManifoldDetection/DensityBasedManifoldDetection.R")
dset <- lapply(paste("./dset", 5:9, ".csv", sep = ""), function(f){
  dset <- as.matrix(read.csv(f))
  rownames(dset) <- as.character(1:nrow(dset))
  dset
})
col <- c(rep(rainbow(8)[1], 100), rep(rainbow(8)[2], 100),
         rep(rainbow(8)[3], 100), rep(rainbow(8)[4], 100),
         rep(rainbow(8)[5], 100), rep(rainbow(8)[6], 100),
         rep(rainbow(8)[7], 100), rep(rainbow(8)[8], 100))
den <- lapply(dset, function(d) sapply(colSums(as.matrix(dist(d)) < 0.5), function(x) max(1, min(x, 20))))
#data preparation

mst <- lapply(dset, function(d) as.matrix(ape::mst(dist(d))))#MST
ddr <- lapply(dset, function(d) as.matrix(DDRTree::DDRTree(t(d), dimensions = 2)$stree) > 0)#Principal Graph
knn <- lapply(dset, function(d){
  lapply(1:12, function(k, d){
    knn <- dbscan::kNN(dist(d), k = k, sort = F)$id
    knn <- cbind(as.numeric(sapply(1:nrow(knn), function(x) rep(x, k))), as.numeric(t(knn)))
    as.matrix(igraph::as_adj(graph_from_edgelist(knn)))
  }, d=d)
})#KNN
dbmd <- lapply(dset, function(d){
  lapply(1:4, function(r, d){
    lapply(c(2, 5, 10), function(n, r, d){
      dbmd <- DensityBasedManifoldDetection(dist(d), radius = r, nobs = n)
      dbmd + t(dbmd) == 2
    }, r=r, d=d)
  }, d=d)
})#DensitiBasedManifoldDetection
#algorithms

pdf("Fig1.pdf", width = 9, height = 6)
layout(matrix(c(1, 1, 1, 2, 4, 6,
                1, 1, 1, 3, 5, 7,
                8, 8, 9, 9, 10, 11,
                8, 8, 9, 9, 12, 12), 4, 6, byrow = T))
par(mar = c(4, 6, 4, 6))
smoothScatter(dset[[5]], nrpoints = 0, xlab = "", ylab = "", axes = F, main = "Simulated Data")
lines(dset[[5]], type = "p", col = col, pch = 16)

par(mar = c(1, 1, 3, 1))
set.seed(1)
plot(igraph::graph_from_adjacency_matrix(mst[[5]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "MST")
set.seed(1)
plot(igraph::graph_from_adjacency_matrix(ddr[[5]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "DDR")
set.seed(1)
plot(igraph::graph_from_adjacency_matrix(knn[[5]][[3]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "KNN, k=3")
set.seed(1)
plot(igraph::graph_from_adjacency_matrix(knn[[5]][[8]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "KNN, k=8")
set.seed(1)
plot(igraph::graph_from_adjacency_matrix(dbmd[[5]][[1]][[3]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "DBMD, r=1, n=10")
set.seed(1)
plot(igraph::graph_from_adjacency_matrix(dbmd[[5]][[2]][[2]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "DBMD, r=2, n=5")

consensus <- mst[[i]] + ddr[[i]] + Reduce("+", knn[[i]])/12 + Reduce("+", lapply(dbmd[[i]], function(x) Reduce("+", x)))/12
svd <- svd(consensus)
D <- diag(svd$d)
diag(D)[50:ncol(D)] <- 0
consensus1 <- svd$u %*% D %*% t(svd$v)
consensus1 <- consensus1 >= 1
set.seed(1)
plot(igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "First 50 SVs")
set.seed(1)
plot(igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected"), vertex.color = topo.colors(20)[den[[5]]], vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "Local Density")

for (j in c(10, 200)){
  D <- diag(svd$d)
  diag(D)[j:ncol(D)] <- 0
  consensus1 <- svd$u %*% D %*% t(svd$v)
  consensus1 <- consensus1 >= 1
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(consensus1, diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = paste("First", j, "SVs"))
}

par(mar = c(3, 5, 1, 1))
plot(svd$d, xlab = "Rank", ylab = "SV", xaxt = "n", yaxt = "n", cex = c(0.5, 2)[1:length(svd$d) %in% c(10, 50, 200)+1], pch = 16, main = "SVD")
axis(side = 1, at = c(10, 50, 200), las = 2)
axis(side = 2, at = seq(0, 60, length.out = 3), las = 2)
dev.off()
#Fig1

pdf("SupFig1.pdf", width = 14, height = 10)
par(mfrow = c(5, 7), mar = c(1, 1, 3, 1))
for (i in 1:length(dset)){
  smoothScatter(dset[[i]], nrpoints = 0, axes = F, xlab = "", ylab = "", main = "Original Dataset")
  lines(dset[[i]], type = "p", col = col, pch = 16, cex = 0.5)
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(mst[[i]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "MST")
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(ddr[[i]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "DDR")
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(knn[[i]][[3]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "KNN, k=3")
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(knn[[i]][[8]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "KNN, k=8")
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(dbmd[[5]][[1]][[3]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "DBMD, r=1, n=10")
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(dbmd[[5]][[2]][[2]], diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "DBMD, r=2, n=5")
}
dev.off()
#SupFig1

pdf("SupFig2.pdf", width = 12, height = 10)
par(mfrow = c(5, 6), mar = c(5, 5, 3, 1))
for (i in 1:length(dset)){
  smoothScatter(dset[[i]], nrpoints = 0, xlab = "", ylab = "", axes = F, main = "Original Dataset")
  lines(dset[[i]], type = "p", col = col, pch = 16, cex = 0.5)
  smoothScatter(den[[i]], degree(igraph::graph_from_adjacency_matrix(mst[[i]], diag = F, mode = "undirected")), xlab = "Local Density", ylab = "Degree", nrpoints = Inf, main = "MST")
  smoothScatter(den[[i]], degree(igraph::graph_from_adjacency_matrix(ddr[[i]], diag = F, mode = "undirected")), xlab = "Local Density", ylab = "Degree", nrpoints = Inf, main = "DDR")
  smoothScatter(den[[i]], degree(igraph::graph_from_adjacency_matrix(knn[[i]][[3]], diag = F, mode = "undirected")), xlab = "Local Density", ylab = "Degree", nrpoints = Inf, main = "KNN, k=3")
  smoothScatter(den[[i]], degree(igraph::graph_from_adjacency_matrix(knn[[i]][[8]], diag = F, mode = "undirected")), xlab = "Local Density", ylab = "Degree", nrpoints = Inf, main = "KNN, k=8")
  smoothScatter(den[[i]], degree(igraph::graph_from_adjacency_matrix(dbmd[[i]][[2]][[2]], diag = F, mode = "undirected")), xlab = "Local Density", ylab = "Degree", nrpoints = Inf, main = "DBMD, r=2, n=5")
}
dev.off()
#SupFig2

pdf("SupFig3.pdf", width = 10, height = 6)
par(mfcol = c(3, 5), mar = c(1, 1, 3, 1))
for (i in 1:length(dset)){
  consensus <- mst[[i]] + ddr[[i]] + Reduce("+", knn[[i]])/12 + Reduce("+", lapply(dbmd[[i]], function(x) Reduce("+", x)))/12
  consensus <- consensus >= 1
  #consensus
  smoothScatter(dset[[i]], nrpoints = 0, xlab = "", ylab = "", axes = F, main = "Original Dataset")
  lines(dset[[i]], type = "p", col = col, pch = 16, cex = 0.5)
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(consensus, diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "Consensus")
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(consensus, diag = F, mode = "undirected"), vertex.color = topo.colors(20)[den[[i]]], vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "Local Density")
}
dev.off()
#SupFig3

overlap <- lapply(1:length(dset), function(j, mst, ddr, knn, dbmd){
  consensus <- mst[[j]] + ddr[[j]] + Reduce("+", knn[[j]])/12 + Reduce("+", lapply(dbmd[[j]], function(x) Reduce("+", x)))/12
  svd <- svd(consensus)
  D <- diag(svd$d)
  diag(D)[50:ncol(D)] <- 0
  consensus <- svd$u %*% D %*% t(svd$v)
  consensus <- consensus >= 1
  
  consensus1 <- mst[[j]] + ddr[[j]] + Reduce("+", knn[[j]])/12
  svd1 <- svd(consensus1)
  D <- diag(svd1$d)
  diag(D)[50:ncol(D)] <- 0
  consensus1 <- svd$u %*% D %*% t(svd$v)
  consensus1 <- consensus1 >= 1

  list(svd=svd, consensus=consensus, svd1=svd1, consensus1=consensus1)
}, mst=mst, ddr=ddr, knn=knn, dbmd=dbmd)
pdf("SupFig4.pdf", width = 10, height = 14)
par(mfcol = c(7, 5))
for (i in 1:length(dset)){
  par(mar = c(1, 1, 3, 1))
  smoothScatter(dset[[i]], nrpoints = 0, xlab = "", ylab = "", axes = F, main = "Original Dataset")
  lines(dset[[i]], type = "p", col = col, pch = 16, cex = 0.5)
  
  par(mar = c(5, 5, 3, 1))
  smoothScatter(den[[i]], degree(igraph::graph_from_adjacency_matrix(overlap[[i]]$consensus, diag = F, mode = "undirected")), xlab = "Local Density", ylab = "Degree", nrpoints = Inf, main = "Consensus")
  smoothScatter(den[[i]], degree(igraph::graph_from_adjacency_matrix(overlap[[i]]$consensus1, diag = F, mode = "undirected")), xlab = "Local Density", ylab = "Degree", nrpoints = Inf, main = "DBMD(-) Consensus")
  
  par(mar = c(1, 1, 3, 1))
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(overlap[[i]]$consensus - mst[[i]] > 0, diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "Consensus - MST")
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(overlap[[i]]$consensus - ddr[[i]] > 0, diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "Consensus - DDR")
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(overlap[[i]]$consensus - (Reduce("+", knn[[i]]) == 12) > 0, diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "Consensus - KNN")
  set.seed(1)
  plot(igraph::graph_from_adjacency_matrix(overlap[[i]]$consensus - (Reduce("+", lapply(dbmd[[i]], function(x) Reduce("+", x))) == 12) > 0, diag = F, mode = "undirected"), vertex.color = col, vertex.frame.color = NA, vertex.label = "", vertex.size = 6, main = "Consensus - DBMD")
}
dev.off()
#SupFig4
