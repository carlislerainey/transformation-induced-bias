
# clear working directory
rm(list = ls())

# set working directory
setwd("~/Dropbox/projects/transformation-induced-bias/")

# load packages
library(compactr)
library(MASS)

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
  x.hi <- qnorm(0.25)
  x.lo <- qnorm(0.75)
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
 mle.me <- mle.pr <- matrix(NA, nrow = n.sims, ncol = length(x0))
  for (i in 1:n.sims) {
    # simulate outcome variable
    y <- rbinom(n, 1, p)
    mle.fit <- glm(y ~ X - 1, family = "binomial")
    # mle
    # mle.fit <- glm(y ~ X - 1, family = "binomial") already done above!
    mle.int[i] <- coef(mle.fit)[1]  
    mle.coef[i] <- coef(mle.fit)[2]
    mle.fd[i] <- plogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x.hi) - 
      plogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x.lo)
    mle.me[i, ] <- dlogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x0)*coef(mle.fit)[2]
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
           true.pr = true.pr,
           e.mle.pr = apply(mle.pr, 2, mean),
           true.rr = true.rr,
           e.mle.rr = mean(mle.rr),
           prop.ones = mean(p))
  return(res)
}

# # test
set.seed(6371)
coef <- int <- fd <- rr <- me <- pr <- 
  total.tau.bias <- coef.tau.bias <- trans.tau.bias <- NULL
n <- round(exp(exp(seq(log(log(100)), log(log(2500)), length.out = 10))))
#n <- c(50, 100, 250,  1000)
for (i in 1:length(n)) {
  cat(paste("\nWorking on N = ", n[i], "...\n", sep = ""))
  sims <- simulate(n = n[i], k = 6, b0 = -1, b1 = 0.15, n.sims = 10000)
  # intercept
  int <- c(int, 100*(sims$e.mle.int - sims$true.int)/sims$true.int)
  # coefficient
  coef <- c(coef, 100*(sims$e.mle.coef - sims$true.coef)/sims$true.coef)
  # total tau bias
  total.tau.bias <- rbind(total.tau.bias, 100*(sims$e.mle.me - sims$true.me)/sims$true.me)
  # coef-induced tau bias
  ctb <- dlogis(sims$e.mle.int + sims$e.mle.coef*sims$x0)*sims$e.mle.coef
  coef.tau.bias <- rbind(coef.tau.bias, 100*(ctb - sims$true.me)/sims$true.me)
  # trans-induced bias
  trans.tau.bias <- rbind(trans.tau.bias, 100*(sims$e.mle.me - ctb)/sims$true.me)
}

pdf("doc/figs/bias-coef.pdf", height = 2.5, width = 6.5)
accept <- 3
xat <- c(100, 250, 500, 1000, 2500, 5000)
ylim <- c(1.12, 1)*mm(c(int, coef))
par(mfrow = c(1, 2), oma = c(3, 3, 1, 1), mar = c(0.5, 0.5, 0.5, 0.5))
eplot(xlim = c(.98, 1.02)*mm(log(n)), ylim = ylim,
      ylab = "Percent Bias",
      xlab = "Sample Size",
      main = expression(hat(beta)[cons]),
      xat = log(xat),
      xticklab = xat
      #yat = c(-10, -5, 0, 5, 10)
      )
# polygon(x = par("usr")[c(1, 2, 2, 1)], 
#         y = c(accept, accept, -accept, -accept),
#         col = "grey90",
#         border = NA)
abline(h = c(accept, -accept), lty = 2, col = "grey50")
abline(h = 0, lty = 1, col = "grey50")
abline(v = log(10*sims$k/sims$prop.ones), lty = 1, col = "grey50")
lines(log(n), int, lwd = 1.7)
#text(log(n[1]), int[1], expression(E(hat(beta)[cons])), pos = 2, cex = 0.6)
text(log(10*sims$k/sims$prop.ones), 
     .9*(par("usr")[4] - par("usr")[3]) + par("usr")[3],
     paste("Rule of thumb: N = ", round(10*sims$k/sims$prop.ones), sep = ""),
     pos = 4, col = "grey50", cex = 0.7)
aplot(expression(hat(beta)[1]))
# polygon(x = par("usr")[c(1, 2, 2, 1)], 
#         y = c(accept, accept, -accept, -accept),
#         col = "grey90",
#         border = NA)
abline(h = c(accept, -accept), lty = 2, col = "grey50")
abline(h = 0, lty = 1, col = "grey50")
abline(v = log(10*sims$k/sims$prop.ones), lty = 1, col = "grey50")
#text()
lines(log(n), coef, lwd = 1.7)
#text(log(n[1]), coef[1], expression(E(hat(beta)[x])), pos = 2, cex = 0.6)
dev.off()

pdf("doc/figs/bias-me.pdf", height = 2.5, width = 9)
par(mfrow = c(1, 3), oma = c(3, 3, 1, 1), mar = c(0.5, 0.5, 0.5, 0.5))
ylim <- c(1.12, 1)*mm(c(total.tau.bias, coef.tau.bias, trans.tau.bias))
eplot(xlim = c(.98, 1.02)*mm(log(n)), ylim = ylim,
      ylab = "Percent Bias",
      xlab = "Sample Size",
      main = expression(paste("Total ", tau, "-bais", sep = "")),
      xat = log(xat),
      xticklab = xat
      #yat = c(-10, -5, 0, 5, 10)
)
# polygon(x = par("usr")[c(1, 2, 2, 1)], 
#         y = c(accept, accept, -accept, -accept),
#         col = "grey90",
#         border = NA)
abline(h = c(accept, -accept), lty = 2, col = "grey50")
abline(h = 0, lty = 1, col = "grey50")
abline(v = log(10*sims$k/sims$prop.ones), lty = 1, col = "grey50")
for (i in 1:ncol(total.tau.bias)) {
  lines(log(n), total.tau.bias[, i], lwd = 1.7)
  text(log(n[1]), total.tau.bias[1, i], bquote(x[1] == .(sims$x0[i])), 
       pos = 2 + sign(total.tau.bias[, i]), 
       cex = 0.6)
}

aplot(expression(paste("Coefficient-Induced ", tau, "-bais", sep = "")))
# polygon(x = par("usr")[c(1, 2, 2, 1)], 
#         y = c(accept, accept, -accept, -accept),
#         col = "grey90",
#         border = NA)
abline(h = c(accept, -accept), lty = 2, col = "grey50")
abline(h = 0, lty = 1, col = "grey50")
abline(v = log(10*sims$k/sims$prop.ones), lty = 1, col = "grey50")
for (i in 1:ncol(total.tau.bias)) {
  lines(log(n), coef.tau.bias[, i], lwd = 1.7)
  #text(log(n[1]), coef.tau.bias[1, i], paste("x = ", sims$x0[i], sep = ""), pos = 2, cex = 0.6)
}

aplot(expression(paste("Transformation-Induced ", tau, "-bais", sep = "")))
# polygon(x = par("usr")[c(1, 2, 2, 1)], 
#         y = c(accept, accept, -accept, -accept),
#         col = "grey90",
#         border = NA)
abline(h = c(accept, -accept), lty = 2, col = "grey50")
abline(h = 0, lty = 1, col = "grey50")
abline(v = log(10*sims$k/sims$prop.ones), lty = 1, col = "grey50")
for (i in 1:ncol(total.tau.bias)) {
  lines(log(n), trans.tau.bias[, i], lwd = 1.7)
  text(log(n[1]), trans.tau.bias[1, i], bquote(x[1] == .(sims$x0[i])), 
       pos = 2 + sign(trans.tau.bias[, i]), 
       cex = 0.6)
}

dev.off()




