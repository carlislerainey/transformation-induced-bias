
# clear working directory
rm(list = ls())

# set working directory
setwd("~/Dropbox/projects/transformation-induced-bias/")

# load packages
library(MASS)
library(scales)
library(ggplot2)

simulate <- function(n, k, b0, b1, n_sims) {
  # set progress bar
  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
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
  # set x_hi and x_lo
  x_hi <- qnorm(0.50)
  x_lo <- qnorm(0.01)
  x0 <- seq(-3, 3, length.out = 20)
  # true values
  true_coef <- b1
  true_int <- b0
  true_fd <- plogis(b0 + b1*x_hi) - plogis(b0 + b1*x_lo)
  true_me <- dlogis(b0 + b1*x0)*b1
  true_pr <- plogis(b0 + b1*x0)
  true_rr <- plogis(b0 + b1*x_hi)/plogis(b0 + b1*x_lo)
  # create holders
 mle_int <- mle_coef <- mle_fd <- mle_rr <- prop_ones <- numeric(n_sims)
 mle_me <- mle_pr <- matrix(NA, nrow = n_sims, ncol = length(x0))
  for (i in 1:n_sims) {
    # simulate outcome variable
    y <- rbinom(n, 1, p)
    mle_fit <- glm(y ~ X - 1, family = "binomial")
    # mle
    # mle_fit <- glm(y ~ X - 1, family = "binomial") already done above!
    mle_int[i] <- coef(mle_fit)[1]  
    mle_coef[i] <- coef(mle_fit)[2]
    mle_fd[i] <- plogis(coef(mle_fit)[1] + coef(mle_fit)[2]*x_hi) - 
      plogis(coef(mle_fit)[1] + coef(mle_fit)[2]*x_lo)
    mle_me[i, ] <- dlogis(coef(mle_fit)[1] + coef(mle_fit)[2]*x0)*coef(mle_fit)[2]
    mle_pr[i, ] <- plogis(coef(mle_fit)[1] + coef(mle_fit)[2]*x0)
    mle_rr[i] <- plogis(coef(mle_fit)[1] + coef(mle_fit)[2]*x_hi)/
      plogis(coef(mle_fit)[1] + coef(mle_fit)[2]*x_lo)
    setTxtProgressBar(pb, i)
  }
  
  res <- list(n = n,
           k = k,
           b0 = b0, 
           b1 = b1, 
           x_hi = x_hi,
           x_lo = x_lo,
           x0 = x0,
           n_sims = n_sims,
           true_int = true_int,
           true_coef = true_coef,
           e_mle_int = mean(mle_int),
           e_mle_coef = mean(mle_coef),
           true_fd = true_fd,
           e_mle_fd = mean(mle_fd),
           true_me = true_me,
           e_mle_me = apply(mle_me, 2, mean),
           true_pr = true_pr,
           e_mle_pr = apply(mle_pr, 2, mean),
           true_rr = true_rr,
           e_mle_rr = mean(mle_rr),
           prop_ones = mean(p))
  return(res)
}

# # test
set.seed(6371)
coef <- int <- fd <- rr <- me <- pr <- total_tau_bias_rr <- 
  total_tau_bias <- coef_tau_bias <- trans_tau_bias <- NULL
n <- round(exp(exp(seq(log(log(100)), log(log(3000)), length.out = 10))))
for (i in 1:length(n)) {
  cat(paste("\nWorking on N = ", n[i], "...\n", sep = ""))
  sims <- simulate(n = n[i], k = 6, b0 = -1, b1 = 0.15, n_sims = 100000)
  # intercept
  int <- c(int, 100*(sims$e_mle_int - sims$true_int)/sims$true_int)
  # coefficient
  coef <- c(coef, 100*(sims$e_mle_coef - sims$true_coef)/sims$true_coef)
  # risk ratio
  total_tau_bias_rr <- c(total_tau_bias_rr, 100*(sims$e_mle_rr - sims$true_rr)/sims$true_rr)

  # total tau bias
  total_tau_bias_df0 <- data.frame(quantity = "'Total'~tau*'-Bias'",
                                   me_where = sims$x0, 
                                   bias = 100*(sims$e_mle_me - sims$true_me)/sims$true_me,
                                   sample_size = n[i])
  total_tau_bias <- rbind(total_tau_bias, total_tau_bias_df0)
  # coef-induced tau bias
  ctb <- dlogis(sims$e_mle_int + sims$e_mle_coef*sims$x0)*sims$e_mle_coef
  coef_tau_bias_df0 <- data.frame(quantity = "'Coefficient-Induced'~tau*'-Bias'",
                                   me_where = sims$x0, 
                                   bias = 100*(ctb - sims$true_me)/sims$true_me,
                                   sample_size = n[i])
  coef_tau_bias <- rbind(coef_tau_bias, coef_tau_bias_df0)
  # trans-induced tau bias
  trans_tau_bias_df0 <- data.frame(quantity = "'Transformation-Induced'~tau*'-Bias'",
                                   me_where = sims$x0, 
                                   bias = 100*(sims$e_mle_me - ctb)/sims$true_me,
                                   sample_size = n[i])
  trans_tau_bias <- rbind(trans_tau_bias, trans_tau_bias_df0)
}

coef_bias_df <- rbind(data.frame(coef = "hat(beta)[cons]", bias = int, sample_size = n), 
                   data.frame(coef = "hat(beta)[1]",  bias = coef, sample_size = n))
me_bias_df <- rbind(total_tau_bias, coef_tau_bias, trans_tau_bias)

gg <- ggplot(coef_bias_df, aes(x = sample_size, y = bias)) + 
  geom_hline(yintercept = c(3, -3), linetype = "dashed") + 
  geom_vline(xintercept = ceiling(60/sims$prop_ones), linetype = "dotted") + 
  annotate(geom = "text", 
           x = ceiling(60/sims$prop_ones), 
           y = Inf, 
           label = paste("N = ", ceiling(60/sims$prop_ones)),
           hjust = -0.1,
           vjust = 1.5,
           size = 3) + 
  geom_line() + 
  facet_grid(~ coef, labeller = "label_parsed") +
  labs(x = "Sample Size",
       y = "Percent Bias",
       title = "Bias in Logistic Regression Coefficients")
ggsave("doc/figs/bias-coef.pdf", gg, height = 3, width = 7)


gg <- ggplot(me_bias_df, aes(x = sample_size, y = bias, group = me_where, color = me_where)) + 
  geom_hline(yintercept = c(3, -3), linetype = "dashed") + 
  geom_vline(xintercept = ceiling(60/sims$prop_ones), linetype = "dotted") + 
  annotate(geom = "text", 
           x = ceiling(60/sims$prop_ones), 
           y = Inf, 
           label = paste("N = ", ceiling(60/sims$prop_ones)),
           hjust = -0.1,
           vjust = 1.5,
           size = 3) + 
  geom_line() + 
  facet_grid(~ quantity, labeller = "label_parsed") +
  scale_color_gradient2(high = "#40004b",
                        mid = "#f7f7f7",
                        low = "#00441b",
                        midpoint = 0) +
  labs(x = "Sample Size",
       y = "Percent Bias",
       title = expression(paste(tau, "-Bias for Marginal Effects")),
       color = expression(x[1]))
ggsave("doc/figs/bias-me.pdf", gg, height = 3, width = 10)


