
# clear workspace
rm(list = ls())

# save simulation data
coef_bias_df <- readRDS("data/mc-sims/logit-sims-coefs.RData")
me_bias_df <- readRDS("data/mc-sims/logit-sims-mes.RData")
sims <- readRDS("data/mc-sims/logit-sims-output.RData")

# plot coefficients
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
       title = "Bias in Logistic Regression Coefficients") + 
  theme_bw()
ggsave("doc/figs/bias-coef.pdf", gg, height = 3, width = 7)

# plot marginal effects
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
  scale_color_gradient2(high = "#f1a340",
                        mid = "#f7f7f7",
                        low = "#998ec3",
                        midpoint = 0) +
  labs(x = "Sample Size",
       y = "Percent Bias",
       title = expression(paste(tau, "-Bias for Marginal Effects")),
       color = expression(x[1])) +
  theme_bw()
ggsave("doc/figs/bias-me.pdf", gg, height = 3, width = 10)


