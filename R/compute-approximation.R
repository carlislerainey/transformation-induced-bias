

# clear working directory
rm(list = ls())

# set seed
set.seed(8742570)

# load packages
library(plot3D)
library(numDeriv)

qi.fn <- function(beta) {
  x0 <- -1
  x1 <- 1
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- plogis(b0 + b1*x1) - plogis(b0 + b1*x0)
  #qi <- dlogis(b0 + b1*x1)*b1
  return(qi)
}

# 3d of beta ll
b0 <- seq(-3, 3, length.out = 50)
b1 <- seq(-3, 3, length.out = 50)
h11 <- h12 <- h21 <- h22 <- ht <- 
  matrix(NA, nrow = length(b0), ncol = length(b1))
for (i in 1:length(b0)) {
  for (j in 1:length(b1)) {
    h <- hessian(qi.fn, c(b0[i], b1[j]))
    h11[i, j] <- h[1, 1]
    h12[i, j] <- h[1, 2]
    h21[i, j] <- h[2, 1]
    h22[i, j] <- h[2, 2]
    true <- qi.fn(c(b0[i], b1[j]))
    ht[i, j] <- 100*sum(diag(h)*c(0.01, 0.01))/true
  }
}

ct <- function(hii) {
  contour(z = hii, x = b0, y = b1,
          xlab = "b0",
          ylab = "b1",
          lwd = 2)
}

par(mfrow = c(3, 3), mar = c(1,1,2,0), oma = c(0,0,0,0))
phi <- c(0, 20, 40)
theta <- c(0, 20, 40)
for (i in 1:length(phi)) {
  for (j in 1:length(theta)) {
    persp(z = ht, x = b0, y = b1,
          xlab = "b0",
          ylab = "b1",
          zlab = "approximate bias",
          phi = phi[i],
          theta = theta[j])
  }
}

#par(mfrow = c(2, 2), mar = c(5,4,1,1), oma = c(0,0,0,0))
#ct(h11); ct(h12); ct(h21); ct(h22)

par(mfrow = c(1, 1), mar = c(5,4,1,1), oma = c(0,0,0,0))
ct(ht)

