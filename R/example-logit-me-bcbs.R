
# clear working directory
rm(list = ls())

# set working directory
setwd("~/Dropbox/projects/transformation-induced-bias/")

# load packages
library(compactr)
library(MASS)
library(boot)

# create function for bootstrap
bs.fn <- function(f, data, indices, x0) {
  data <- data[indices,] # allows boot to select sample 
  bs.fit <- glm(f, data, family = "binomial")
  me <- dlogis(coef(bs.fit)[1] + coef(bs.fit)[2]*x0)*coef(bs.fit)[2]
  return(me) 
} 

simulate <- function(n, k, b0, b1, n.sims) {
  # set progress bar
  pb <- txtProgressBar(min = 0, max = n.sims, style = 3)
  # set covariates and probabilty
  rho <- 0.0
  Sigma <- matrix(rho, nrow = k, ncol = k)
  diag(Sigma) <- 1
  X <- mvrnorm(n, mu = rep(0, k), Sigma = Sigma)
  X <- cbind(1, X)
  #x <- round(seq(10, 20, length.out = n))
  #X <- cbind(1, x)
  b <- c(b0, b1, rep(.15, k - 1))
  p <- plogis(X%*%b)
  # set x.hi and x.lo
  x.hi <- qnorm(1.00)
  x.lo <- qnorm(0.50)
  x0 <- c(-3, -2, -1, 0, 1, 2, 3)#seq(-3, 3, length.out = 20)
  # true values
  true.coef <- b1
  true.int <- b0
  true.fd <- plogis(b0 + b1*x.hi) - plogis(b0 + b1*x.lo)
  true.me <- dlogis(b0 + b1*x0)*b1
  true.pr <- plogis(b0 + b1*x0)
  true.rr <- plogis(b0 + b1*x.hi)/plogis(b0 + b1*x.lo)
  # create holders
 mle.int <- mle.coef <- mle.fd <- mle.rr <- prop.ones <- numeric(n.sims)
 mle.me <- mle.me.bc <- mle.pr <- matrix(NA, nrow = n.sims, ncol = length(x0))
  for (i in 1:n.sims) {
    # simulate outcome variable
    y <- rbinom(n, 1, p)
    d <- data.frame(y, X)
    mle.fit <- glm(y ~ X - 1, family = "binomial")
    bs.fit <- boot(data = d, bs.fn, R = 100, ncpus = 4,
                   x0 = x0, f = y ~ X2 + X3 + X4 + X5)
    bias <- apply(bs.fit$t, 2, mean) - bs.fit$t0
    # mle
    # mle.fit <- glm(y ~ X - 1, family = "binomial") already done above!
    mle.int[i] <- coef(mle.fit)[1]  
    mle.coef[i] <- coef(mle.fit)[2]
    mle.fd[i] <- plogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x.hi) - 
      plogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x.lo)
    mle.me[i, ] <- dlogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x0)*coef(mle.fit)[2]
    mle.me.bc[i, ] <- dlogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x0)*coef(mle.fit)[2] - bias
    mle.pr[i, ] <- plogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x0)
    mle.rr[i] <- plogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x.hi)/
      plogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x.lo)
    setTxtProgressBar(pb, i)
  }
  
  res <- list(n = n,
           k = k,
           b0 = b0, 
           b1 = b1, 
           x.hi = x.hi,
           x.lo = x.lo,
           x0 = x0,
           n.sims = n.sims,
           true.int = true.int,
           true.coef = true.coef,
           e.mle.int = mean(mle.int),
           e.mle.coef = mean(mle.coef),
           true.fd = true.fd,
           e.mle.fd = mean(mle.fd),
           true.me = true.me,
           e.mle.me = apply(mle.me, 2, mean),
           e.mle.me.bc = apply(mle.me.bc, 2, mean),
           true.pr = true.pr,
           e.mle.pr = apply(mle.pr, 2, mean),
           true.rr = true.rr,
           e.mle.rr = mean(mle.rr),
           prop.ones = mean(p))
  return(res)
}


sims <- simulate(n = 500, k = 4, b0 = -1, b1 = 0.15, n.sims = 1000)

true <- sims$true.me
me <- sims$e.mle.me
bc <- sims$e.mle.me.bc

round(100*(me - true)/true)
round(100*(bc - true)/true)


