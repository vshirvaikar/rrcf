# Relative Risk Causal Forests

A modification of the <a href="https://github.com/grf-labs/grf">Generalized Random Forests (grf)</a>
package that targets relative risk heterogeneity in causal forests. The **simulations** folder contains
all code necessary to replicate the experiments in Shirvaikar, Storås, Lin and Holmes (2024).

### Installation

This development code can be installed from source using devtools. 
The original GRF package should be installed first, as it provides all analysis 
functions related to causal forests (variable importance, visualization, etc.)

```R
install.packages("grf")
devtools::install_github("vshirvaikar/rrcf", subdir = "r-package/rrcf")

library(grf)
library(rrcf)
```

***Note***: As of 2025, the installation is supported for MacOS and Linux, but runs into a
compilation error on Windows. We are aware of the issue and looking into a solution.

### Guidance

The following script demonstrates how to implement a relative risk causal forest.

```R
set.seed(123)
n <- 2000
p <- 10
X <- as.data.frame(matrix(rnorm(n * p), n, p))
names(X) <- paste0("X", 1:p)

# Generate treatment (with optional confounding)
rct_flag <- TRUE # set to FALSE for observational data
if (rct_flag) {
  W <- rbinom(n, 1, 0.5)
} else {
  W <- rbinom(n, 1, plogis(X$X1 - X$X2))
}

# Generate binary outcome with heterogeneous relative risk
p0 <- plogis(-1 + 0.6 * X$X2)       # baseline risk
log_rr <- pmax(X$X1, 0)             # heterogeneity in log-RR
p1 <- pmin(p0 * exp(log_rr), 0.95)  # treated risk = RR * baseline (capped)
Y <- rbinom(n, 1, ifelse(W == 1, p1, p0))

# Train/test split
idx <- sample.int(n, floor(0.8 * n))
X.train <- X[idx, , drop = FALSE]
Y.train <- Y[idx]
W.train <- W[idx]
X.test  <- X[-idx, , drop = FALSE]
Y.test  <- Y[-idx]
W.test  <- W[-idx]

# Fit relative risk causal forest
forest.rr <- rr_causal_forest(X.train, Y.train, W.train, rct = rct_flag)

# Predict relative conditional treatment effect
rr.pred <- rr_predict(forest.rr, X.test)

# Test for treatment effect heterogeneity
pval <- rr_test_calibration(forest.rr, X.test, Y.test, W.test)
```

### References

Vik Shirvaikar, Andrea Storås, Xi Lin and Chris Holmes.
<b>Targeting relative risk heterogeneity with causal forests.</b> <i>Submitted</i>, 2024. 
[<a href="https://arxiv.org/abs/2309.15793">arXiv</a>]

Susan Athey, Julie Tibshirani and Stefan Wager.
<b>Generalized Random Forests.</b> <i>Annals of Statistics</i>, 47(2), 2019.
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>,
<a href="https://arxiv.org/abs/1610.01271">arXiv</a>]

Stefan Wager and Susan Athey.
<b>Estimation and Inference of Heterogeneous Treatment Effects using Random Forests.</b>
<i>Journal of the American Statistical Association</i>, 113(523), 2018.
[<a href="https://www.tandfonline.com/eprint/v7p66PsDhHCYiPafTJwC/full">paper</a>,
<a href="https://arxiv.org/abs/1510.04342">arXiv</a>]
