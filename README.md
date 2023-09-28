# Relative Risk Causal Forests

A modification of the <a href="https://github.com/grf-labs/grf">Generalized Random Forests (grf)</a>
package that targets relative risk heterogeneity in causal forests.

### Installation

This development code can be installed from source using devtools. The original GRF package should be installed first, as it provides all analysis functions related to causal forests (variable importance, visualization, etc.)

```R
install.packages("grf")
devtools::install_github("vshirvaikar/rrcf", subdir = "r-package/rrcf")

library(grf)
library(rrcf)
```

### Guidance

The following script demonstrates how to implement a relative risk causal forest.

```R
# Generate data.
n <- 2000
p <- 10
X <- matrix(rnorm(n * p), n, p)
X.test <- matrix(0, 101, p)
X.test[, 1] <- seq(-2, 2, length.out = 101)
W <- rbinom(n, 1, 0.4 + 0.2 * (X[, 1] > 0))
Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

# Train a relative risk causal forest.
forest.rr = rr_causal_forest(X, Y, W)

# View predictions of RELATIVE conditional treatment effect.
predict.rr = rr_predict(forest.rr, X.test)
```

### References

Vik Shirvaikar and Chris Holmes.
<b>Targeting Relative Risk Heterogeneity with Causal Forests.</b> <i>Submitted</i>, 2023. 
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