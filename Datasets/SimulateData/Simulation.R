library(MASS)
library(scatterplot3d)

x <- seq(0, pi, length.out = 100)
y <- sin(x)
dset1 <- cbind(x, y)
write.csv(dset1, file = "dset1.csv", quote = F, row.names = F)
#linear sin

x <- sin(seq(0, 2*pi, length.out = 100))
y <- cos(seq(0, 2*pi, length.out = 100))
z <- sin(seq(0, 4*pi, length.out = 100))
dset2 <- cbind(x, y, z)
write.csv(dset2, file = "dset2.csv", quote = F, row.names = F)
#twisted circle

x <- c(sin(seq(0, 0.5*pi, length.out = 50)), 2 - sin(seq(0, 0.5*pi, length.out = 50)), rep(1, 100))
y <- c(cos(seq(0, 0.5*pi, length.out = 50)), cos(seq(0, 0.5*pi, length.out = 50)), seq(-1, 0, length.out = 100))
dset3 <- cbind(x, y)
write.csv(dset3, file = "dset3.csv", quote = F, row.names = F)
#bifurcation

x <- c(sin(seq(0, 0.5*pi, length.out = 50)),
       sin(seq(0, 0.5*pi, length.out = 50)),
       2 - sin(seq(0, 0.5*pi, length.out = 50)),
       2 - sin(seq(0, 0.5*pi, length.out = 50)),
       rep(1, 50))
y <- c(cos(seq(0, 0.5*pi, length.out = 50)),
       1 + cos(seq(0, 0.5*pi, length.out = 50)),
       cos(seq(0, 0.5*pi, length.out = 50)),
       1 + cos(seq(0, 0.5*pi, length.out = 50)),
       seq(-1, 1, length.out = 50))
dset4 <- cbind(x, y)
write.csv(dset4, file = "dset4.csv", quote = F, row.names = F)
#tree

c1 <- mvrnorm(n = 100, c(-5, -3), matrix(c(5, 3, 3, 2), 2, 2))
c1 <- c1[order(c1[, 2]), ]
c2 <- mvrnorm(n = 100, c(5, 3), matrix(c(5, 3, 3, 2), 2, 2))
c2 <- c2[order(c2[, 2]), ]
dset5 <- rbind(c1, c2)
colnames(dset5) <- c("x", "y")
write.csv(dset5, file = "dset5.csv", quote = F, row.names = F)
#density-linear

c1 <- mvrnorm(n = 100, c(5, 0), matrix(c(0.1, 0, 0, 8), 2, 2))
c1 <- c1[order(c1[, 2]), ]
c2 <- mvrnorm(n = 100, c(0, -3), matrix(c(5, -3, -3, 2), 2, 2))
c2 <- c2[order(c2[, 2]), ]
c3 <- mvrnorm(n = 100, c(0, 3), matrix(c(5, 3, 3, 2), 2, 2))
c3 <- c3[order(c3[, 2]), ]
dset6 <- rbind(c1, c2, c3)
colnames(dset6) <- c("x", "y")
write.csv(dset6, file = "dset6.csv", quote = F, row.names = F)
#density-circle

c1 <- mvrnorm(n = 100, c(0, 12), matrix(c(0.1, 0, 0, 10), 2, 2))
c1 <- c1[order(c1[, 2]), ]
c2 <- mvrnorm(n = 100, c(0, 0), matrix(c(0.1, 0, 0, 10), 2, 2))
c2 <- c2[order(c2[, 2]), ]
c3 <- mvrnorm(n = 100, c(-5, -10), matrix(c(5, 3, 3, 2), 2, 2))
c3 <- c3[order(c3[, 2]), ]
c4 <- mvrnorm(n = 100, c(5, -10), matrix(c(5, -3, -3, 2), 2, 2))
c4 <- c4[order(c4[, 2]), ]
dset7 <- rbind(c1, c2, c3, c4)
colnames(dset7) <- c("x", "y")
write.csv(dset7, file = "dset7.csv", quote = F, row.names = F)
#density-bifurcation

c1 <- mvrnorm(n = 100, c(0, 0), matrix(c(0.1, 0, 0, 10), 2, 2))
c1 <- c1[order(c1[, 2]), ]
c2 <- mvrnorm(n = 100, c(-5, -10), matrix(c(5, 3, 3, 2), 2, 2))
c2 <- c2[order(c2[, 2]), ]
c3 <- mvrnorm(n = 100, c(5, -10), matrix(c(5, -3, -3, 2), 2, 2))
c3 <- c3[order(c3[, 2]), ]
c4 <- mvrnorm(n = 100, c(-5, 10), matrix(c(5, -3, -3, 2), 2, 2))
c4 <- c4[order(c4[, 2]), ]
c5 <- mvrnorm(n = 100, c(5, 10), matrix(c(5, 3, 3, 2), 2, 2))
c5 <- c5[order(c5[, 2]), ]
dset8 <- rbind(c1, c2, c3, c4, c5)
colnames(dset8) <- c("x", "y")
write.csv(dset8, file = "dset8.csv", quote = F, row.names = F)
#density-tree

c1 <- cbind(5*(sin(seq(0, 2*pi, length.out = 100)) + rnorm(100, 0, 0.1) + 10),
            5*cos(seq(0, 2*pi, length.out = 100)) + rnorm(100, 0, 0.1))
c1 <- c1[order(c1[, 2]), ]
c2 <- mvrnorm(n = 100, c(0, 0), matrix(c(0.1, 0, 0, 10), 2, 2))
c2 <- c2[order(c2[, 2]), ]
c3 <- mvrnorm(n = 100, c(-5, -10), matrix(c(5, 3, 3, 2), 2, 2))
c3 <- c3[order(c3[, 2]), ]
c4 <- mvrnorm(n = 100, c(5, -10), matrix(c(5, -3, -3, 2), 2, 2))
c4 <- c4[order(c4[, 2]), ]
c5 <- mvrnorm(n = 100, c(-5, 10), matrix(c(5, -3, -3, 2), 2, 2))
c5 <- c5[order(c5[, 2]), ]
c6 <- mvrnorm(n = 100, c(5, 10), matrix(c(5, 3, 3, 2), 2, 2))
c6 <- c6[order(c6[, 2]), ]
dset9 <- rbind(c1, c2, c3, c4, c5, c6)
colnames(dset9) <- c("x", "y")
write.csv(dset9, file = "dset9.csv", quote = F, row.names = F)
#density-tree-cluster

pdf("SimulatedData.pdf", width = 10, height = 10)
par(mfrow = c(3, 3))
plot(dset1, xlab = "x", ylab = "y", pch = 16, main = "dset1, Linear")
scatterplot3d(dset2, xlab = "x", ylab = "y", zlab = "z", pch = 16, main = "dset2, Circle")
plot(dset3, xlab = "x", ylab = "y", pch = 16, main = "dset3, Bifurcation")
plot(dset4, xlab = "x", ylab = "y", pch = 16, main = "dset4, Tree")
smoothScatter(dset5, xlab = "x", ylab = "y", nrpoints = Inf, main = "dset5, Density-Linear")
smoothScatter(dset6, xlab = "x", ylab = "y", nrpoints = Inf, main = "dset6, Density-Circle")
smoothScatter(dset7, xlab = "x", ylab = "y", nrpoints = Inf, main = "dset7, Density-Bifurcation")
smoothScatter(dset8, xlab = "x", ylab = "y", nrpoints = Inf, main = "dset8, Density-Tree")
smoothScatter(dset9, xlab = "x", ylab = "y", nrpoints = Inf, main = "dset9, Density-Tree-Cluster")
dev.off()
