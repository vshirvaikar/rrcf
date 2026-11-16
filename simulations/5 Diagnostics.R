# Plot power curves - export 6.5x3in
pvals = data.frame(pvals)
colnames(pvals) = c("n", "rho", "GRF", "RRCF")

power = pvals %>% group_by(n, rho) %>%
  summarise(across(c(GRF, RRCF), ~ mean(. < 0.05))) %>%
  pivot_longer(cols = c(GRF, RRCF), names_to = "Method", values_to = "power")

ggplot(power, aes(x = n, y = power, color = factor(rho), linetype = Method)) +
  geom_line(size = 1) + theme_minimal() + theme(legend.position = "right") +
  scale_linetype_manual(values = c("GRF" = "dotted", "RRCF" = "solid")) +
  labs(x = "Sample Size (n)", y = "Power (% of Simulations)", 
       color = "Heterogeneity (rho)")

# Plot error curves - export 6.5x3in
mapes = data.frame(mapes)
colnames(mapes) = c("n", "rho", "GRF", "RRCF")

error = mapes %>% group_by(n, rho) %>%
  summarise(across(c(GRF, RRCF), mean)) %>%
  pivot_longer(cols = c(GRF, RRCF), names_to = "Method", values_to = "mape")

ggplot(error, aes(x = n, y = mape, color = factor(rho), linetype = Method)) +
  geom_line(size = 1) + theme_minimal() + theme(legend.position = "right") +
  scale_linetype_manual(values = c("GRF" = "dotted", "RRCF" = "solid")) +
  labs(x = "Sample Size (n)", y = "Mean Absolute Percentage Error", 
       color = "Heterogeneity (rho)")

# Plot variable importance curves - export 6.5x3in
vimps = data.frame(vimps)
colnames(vimps) = c("trial", "n", "rho", "feature_idx", "GRF", "RRCF")

prog = vimps %>% filter(feature_idx %in% 4:7) %>% group_by(trial, n, rho) %>%
  summarise(across(c(GRF, RRCF), sum), .groups = "drop") %>% group_by(n, rho) %>%
  summarise(across(c(GRF, RRCF), mean), .groups = "drop") %>%
  pivot_longer(cols = c(GRF, RRCF), names_to = "Method", values_to = "vimp")

ggplot(prog, aes(x = n, y = vimp, color = factor(rho), linetype = Method)) +
  geom_line(size = 1) + theme_minimal() + theme(legend.position = "right") +
  scale_linetype_manual(values = c("GRF" = "dotted", "RRCF" = "solid")) +
  labs(x = "Sample Size (n)", y = "Importance of Predictive Variables", 
       color = "Heterogeneity (rho)")

# Summary metrics and analysis for tables
pvals_nn = subset(pvals, rho > 0)
mapes_nn = subset(mapes, rho > 0)
prog_nn = subset(prog, rho > 0)

100 * (mean(pvals_nn$RRCF < 0.05) - mean(pvals_nn$GRF < 0.05))
mean(mapes_nn$GRF) - mean(mapes_nn$RRCF)
mean(subset(prog_nn, Method == "RRCF")$vimp) - mean(subset(prog_nn, Method == "GRF")$vimp)

results_table = power %>% select(n, rho, Method, power) %>%
  left_join(error %>% select(n, rho, Method, mape), by = c("n", "rho", "Method")) %>%
  left_join(prog %>% select(n, rho, Method, vimp), by = c("n", "rho", "Method")) %>%
  pivot_wider(names_from = Method, values_from = c(power, mape, vimp)) %>% arrange(rho, n)
results_table

# View splits per tree to diagnose forest issues
splits_per_tree <- function(forest){return(sapply(1:forest[['_num_trees']], 
                   function(i) length(forest[['_split_vars']][[i]])))}
table(splits_per_tree(forest.grf))
table(splits_per_tree(forest.glm))

