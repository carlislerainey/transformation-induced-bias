

# clear working directory
rm(list = ls())

# set seed
set.seed(8742570)

# load packages
library(numDeriv)
library(dplyr)
library(ggplot2)

# predicted probabilities
qi_pr_lo <- function(beta) {
  x1 <- -2
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- plogis(b0 + b1*x1)  # probability
  return(qi)
}
qi_pr_mid <- function(beta) {
  x1 <- 0
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- plogis(b0 + b1*x1)  # probability
  return(qi)
}
qi_pr_hi <- function(beta) {
  x1 <- 2
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- plogis(b0 + b1*x1)  # probability
  return(qi)
}
# marginal effects
qi_me_lo <- function(beta) {
  x1 <- -2
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- dlogis(b0 + b1*x1)*b1 # marginal effect
  return(qi)
}
qi_me_mid <- function(beta) {
  x1 <- 0
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- dlogis(b0 + b1*x1)*b1 # marginal effect
  return(qi)
}
qi_me_hi <- function(beta) {
  x1 <- 2
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- dlogis(b0 + b1*x1)*b1 # marginal effect
  return(qi)
}
# risk ratio
qi_rr_mod <- function(beta) {
  x0 <- -1
  x1 <- 1
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- plogis(b0 + b1*x1)/plogis(b0 + b1*x0)  # risk ratio
  return(qi)
}
qi_rr_extreme <- function(beta) {
  x0 <- -2
  x1 <- 2
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- plogis(b0 + b1*x1)/plogis(b0 + b1*x0)  # risk ratio
  return(qi)
}
# first difference
qi_fd_mod <- function(beta) {
  x0 <- -1
  x1 <- 1
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- plogis(b0 + b1*x1) - plogis(b0 + b1*x0)  # first difference
  return(qi)
}
qi_fd_extreme <- function(beta) {
  x0 <- -1
  x1 <- 1
  b0 <- beta[1] 
  b1 <- beta[2]
  qi <- plogis(b0 + b1*x1) - plogis(b0 + b1*x0)  # first difference
  return(qi)
}

qi_fn <- list(qi_pr_lo, qi_pr_mid, qi_pr_hi, qi_me_lo, qi_me_mid, qi_me_hi,
              qi_rr_mod, qi_rr_extreme, 
              qi_fd_mod, qi_fd_extreme)
qi_name <- c("Predicted Probability for Low Value of X",
             "Predicted Probability for Middle Value of X",
             "Predicted Probability for High Value of X",
             "Marginal Effect for Low Value of X",
             "Marginal Effect for Middle Value of X",
             "Marginal Effect for High Value of X",
             "Risk Ratio for Moderate Values of X",
             "Risk Ratio for Extreme Values of X",
             "First Difference for Moderate Values of X",
             "First Difference for Extreme Values of X")
qi_type <- c("pp", "pp", "pp", 
             "me", "me", "me", 
             "rr", "rr", 
             "fd", "fd")


b0 <- seq(-3, 3, length.out = 20)
b1 <- seq(-3, 3, length.out = 20)
cov_wt <- c(0, .25, .5)
df <- data.frame(expand.grid(b0, b1, cov_wt, qi_name))
names(df) <- c("b0", "b1", "cov_wt", "qi_name")
df <- left_join(df, data.frame(qi_name = qi_name, qi_type = qi_type))
df$qi_name <- factor(df$qi_name, levels = qi_name)

for (i in 1:length(b0)) {
  for (j in 1:length(b1)) {
    for (q in 1:length(qi_fn)) {
      h_qi <- hessian(qi_fn[[q]], c(b0[i], b1[j]))
      for (k in 1:length(cov_wt)) {
        sum <- sum(diag(h_qi)) + cov_wt[k]*(sum(h_qi) - sum(diag(h_qi)))
        df$sum_hes[df$b0 == b0[i] & 
                     df$b1 == b1[j] & 
                     df$cov_wt == cov_wt[k] &
                     df$qi_name == qi_name[q]] <- sum
      }
    }
  }
}

df <- subset(df, cov_wt == 0)

gg <- ggplot(subset(df, qi_type == "pp"), aes(x = b0, y = b1, z = sum_hes, color = ..level..)) + 
  #geom_tile(aes(fill = sum_hes), alpha = 0.5) + stat_contour() +
  stat_contour(bins = 50) +
  #stat_contour(bins = 10, size = 1.4) +
  facet_wrap(~ qi_name) + 
  scale_color_gradient2(high = "#40004b",
                        mid = "#f7f7f7",
                        low = "#00441b",
                        midpoint = 0) +
  labs(x = expression(beta[0]),
       y = expression(beta[1]),
       color = "Total Curvature",
       title = "Total Diagonal Curvature in the Transformation to Predicted Pobabilities")
ggsave("doc/figs/trans-pr.pdf", gg, height = 5, width = 12)

gg <- gg %+% subset(df, qi_type == "me") + 
  labs(title = "Total Diagonal Curvature in the Transformation to Marginal Effects")
ggsave("doc/figs/trans-me.pdf", gg, height = 5, width = 12)

gg <-gg %+% subset(df, qi_name == "Risk Ratio for Moderate Values of X") + 
  labs(title = "Total Diagonal Curvature in the Transformation to Moderate Risk Ratios")
ggsave("doc/figs/trans-rr1.pdf", gg, height = 4, width = 7)
gg <- gg %+% subset(df, qi_name == "Risk Ratio for Extreme Values of X") + 
  labs(title = "Total Diagonal Curvature in the Transformation to Extreme Risk Ratios")
ggsave("doc/figs/trans-rr2.pdf", gg, height = 4, width = 7)

gg <- gg %+% subset(df, qi_type == "fd") + 
  labs(title = "Total Diagonal Curvature in the Transformation to First Differences")
ggsave("doc/figs/trans-fd.pdf", gg, height = 4, width = 8.5)




