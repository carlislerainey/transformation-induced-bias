
# clear workspace
rm(list = ls())

# read data
df <- readRDS("data/mc-sims/lacina-sims.RData")

# expected value
gg1 <- ggplot(filter(df, qi == "Expected~Value"), aes(x = true, y = bias/true)) + 
  geom_point() + 
  #facet_grid(~ type, labeller = "label_parsed") + 
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent) + 
  labs(title = "Expected Value",
       x = "True Value",
       y = "Percent Bias") + 
  theme_bw()
  

# first difference
gg2 <- gg1 %+% subset(df, qi == "First~Difference") + 
  labs(title = "First Difference")

# combine plots
pdf("doc/figs/lacina.pdf", width = 9, height = 3)
grid.arrange(gg1, gg2, ncol = 2)
dev.off()
