library(tidyr)

# View nice plot of p-values
ggplot(data.frame(pvals05), aes(x = -log10(X2), y = -log10(X1))) +
  theme_minimal() + geom_point() + geom_abline(intercept = 0, slope = 1) +
  geom_vline(xintercept = -log10(0.1), color = "red") + #xlim(0, 2.5) + 
  geom_hline(yintercept = -log10(0.1), color = "red") + #ylim(0, 2.5) + 
  labs(x = "Relative risk forest", y = "GRF causal forest")
colSums(p.values < 0.01)

splits_per_tree <- function(forest){return(sapply(1:forest[['_num_trees']], 
                   function(i) length(forest[['_split_vars']][[i]])))}
table(splits_per_tree(forest.grf))
table(splits_per_tree(forest.glm))

# Get summary metrics
(colSums(pvals < 0.05)[4]-colSums(pvals < 0.05)[3])/12
(colMeans(mapes)[4]-colMeans(mapes)[3])/colMeans(mapes)[3]

# Plot power curves - export 6.5x3in and crop
power = data.frame(pvals) %>% group_by(X1, X2) %>% 
  summarise(across(c(X3, X4), ~ sum(. < 0.05))) %>%   
  rename(GRF = X3, RRCF = X4) %>%
  pivot_longer(cols = c(GRF, RRCF), names_to = "Method", values_to = "power")

ggplot(power, aes(x = X1, y = power, color = factor(X2), linetype = Method)) +
  geom_line(size = 1) + theme_minimal() + theme(legend.position = "right") +
  scale_x_continuous(breaks = c(2500, 5000, 7500, 10000)) +
  scale_linetype_manual(values = c("GRF" = "dotted", "RRCF" = "solid")) +
  labs(x = "Sample Size (n)", y = "Power (%)", color = "Heterogeneity (rho)")

# Plot error curves - export 6.5x3in and crop
error = data.frame(mapes) %>% group_by(X1, X2) %>% 
  summarise(GRF = mean(X3), RRCF = mean(X4)) %>%
  pivot_longer(cols = c(GRF, RRCF), names_to = "Method", values_to = "error")

ggplot(error, aes(x = X1, y = error, color = factor(X2), linetype = Method)) +
  geom_line(size = 1) + theme_minimal() + theme(legend.position = "right") +
  scale_x_continuous(breaks = c(2500, 5000, 7500, 10000)) +
  scale_linetype_manual(values = c("GRF" = "dotted", "RRCF" = "solid")) +
  labs(x = "Sample Size (n)", y = "Average MAPE", color = "Heterogeneity (rho)")
