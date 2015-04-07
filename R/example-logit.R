
# clear workspace
rm(list = ls())

# load packages
library(compactr)
library(logistf)
library(MASS)

# create function to do the simulations
simulate <- function(n, b0, b1, n.sims) {
  # set progress bar
  pb <- txtProgressBar(min = 0, max = n.sims, style = 3)
  # set covariates and probabilty
  b <- c(b0, b1)
  x <- round(seq(10, 20, length.out = n)) #round(runif(n, 10, 21))
  X <- cbind(1, x)
  p <- plogis(X%*%b)
  n.me <- 200
  x0 <- sort(unique(x))
  # true values
  true.me <- dlogis(b0 + b1*x0)*b1
  # create holders
  lpm.me <- mle.me <- pmle.me <- matrix(NA, nrow = n.sims, ncol = length(x0))
  mle.coef <- pmle.coef <- matrix(NA, nrow = n.sims, ncol = 2)

  for (i in 1:n.sims) {
    # simulate outcome variable
    y <- rbinom(n, 1, p)
    mle.fit <- glm(y ~ x, family = "binomial")
    mle.coef[i, ] <- coef(mle.fit)
    # just start over if data set has separation
    #if (sum(abs(coef(mle.fit)) > 7) == 0) {
      # lpm
      lpm.fit <- lm(y ~ x)
      lpm.me[i, ] <- coef(lpm.fit)[2]
      # mle
      mle.me[i,] <- dlogis(coef(mle.fit)[1] + coef(mle.fit)[2]*x0)*coef(mle.fit)[2]
      # pmle
      pmle.fit <- logistf(y ~ x)
      pmle.coef[i, ] <- coef(pmle.fit)
      pmle.me[i, ] <- dlogis(coef(pmle.fit)[1] + coef(pmle.fit)[2]*x0)*coef(pmle.fit)[2]
    #}
    setTxtProgressBar(pb, i)
  }
  res <- list(n = n,
              x = x, 
              b0 = b0, 
              b1 = b1, 
              x0 = x0, 
              p = p,
              n.sims = n.sims,
              true.me = true.me,
              lpm.me = lpm.me,
              mle.me = mle.me,
              pmle.me = pmle.me,
              mle.coef = mle.coef,
              pmle.coef = pmle.coef)
  return(res)
}

# do the simulation
set.seed(4367)
sims <- simulate(n = 30, b0 = -2.5, b1 = .2, n.sims = 100000)

par(mfrow = c(1, 1))
plot(sims$x, sims$p)

# plot the average mes against the true mes
pdf("doc/figs/logit-me-bias.pdf", height = 4, width = 6)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0),
    mar = c(3, 4, 1, 1))
trans <- .8
ylim <- range(c(apply(sims$lpm.me, 2, mean, na.rm = TRUE),
                apply(sims$mle.me, 2, mean, na.rm = TRUE),
                apply(sims$pmle.me, 2, mean, na.rm = TRUE)))
eplot(xlim = mm(sims$x0), ylim = ylim, 
      xlab = "Education", 
      ylab = "Marginal Effect of Education on Pr(y = 1)",
      ylabpos = 2.8)
e.mle.me <- apply(sims$mle.me, 2, mean, na.rm = TRUE)
e.pmle.me <- apply(sims$pmle.me, 2, mean, na.rm = TRUE)
# for (i in 1:length(sims$x0)) {
#   true.me <- sims$true.me[i]
#   red.me <- e.mle.me[i]
#   blue.me <- e.pmle.me[i]
#   which.min <- which.min(c(abs(red.me - true.me), 
#                            abs(blue.me - true.me)))
#   if (which.min == 1) {
#     lines(c(sims$x0[i], sims$x0[i]),
#           c(true.me, red.me), 
#           col = rgb(1, 0, 0, trans))
#   }
#   if (which.min == 2) {
#     lines(c(sims$x0[i], sims$x0[i]),
#           c(true.me, blue.me), 
#           col = rgb(0, 0, 1, trans))
#   }
# }
points(sims$x0, sims$true.me,
     ylim = ylim, pch = 19)
points(sims$x0, e.mle.me,
       ylim = ylim, pch = 1, col = rgb(1, 0, 0, trans))
points(sims$x0, e.pmle.me,
       ylim = ylim, pch = 2, col = rgb(0, 0, 1, trans))
# points(sims$x0, apply(sims$lpm.me, 2, mean, na.rm = TRUE),
#        ylim = ylim, col = 4)

text(10, sims$true.me[1], labels = "true", pos = 3, cex = .7)
text(10, e.mle.me[1], labels = "mle", pos = 3, cex = .7, col = "red")
text(10, e.pmle.me[1], labels = "pmle", pos = 1, cex = .7, col = "blue")
dev.off()




# calculate the bias
apply(sims$mle.coef, 2, mean)/c(-2.5, 0.2)
apply(sims$pmle.coef, 2, mean)/c(-2.5, 0.2)
apply(sims$mle.me, 2, mean)
apply(sims$pmle.me, 2, mean)
sims$true.me

# ratio of average marginal effects
a1 <- sum(table(sort(sims$x))*sims$true.me)/sims$n
a2 <- sum(table(sort(sims$x))*e.mle.me)/sims$n
a3 <- sum(table(sort(sims$x))*e.pmle.me)/sims$n
a1/a2
a1/a3




