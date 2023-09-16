# Relative Risk Causal Forests

PRELIMINARY RELEASE

A modification of the GRF package that targets relative risk heterogeneity in causal forests.

### Installation

This development code can be installed from source using devtools.

```R
devtools::install_github("vshirvaikar/rrcf", subdir = "r-package/rrcf")
```

### Guidance

The following script demonstrates how to implement a relative risk causal forest. As this code is preliminary, it uses some workarounds that are still in the process of being streamlined.

```R
library(rrcf)

# Generate data.
n <- 2000
p <- 10
X <- matrix(rnorm(n * p), n, p)
X.test <- matrix(0, 101, p)
X.test[, 1] <- seq(-2, 2, length.out = 101)
W <- rbinom(n, 1, 0.4 + 0.2 * (X[, 1] > 0))
Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

# Train a relative risk causal forest.
# 1. Since the forest splitting rule is based on Poisson regression, the Y and W
#    inputs should not be centered as in the original GRF code, so the Y.hat and
#    W.hat inputs are instead specified as vectors of zeros.
# 2. The imbalance.penalty value of 100 serves as a flag within the code to use
#    the relative risk splitting rule. This is in the process of being spun out
#    as a separate rr_forest function.
zeros = numeric(length(Y))
forest.rr = causal_forest(X, Y, W, Y.hat=zeros, W.hat=zeros, imbalance.penalty=100)
```

### References

Vik Shirvaikar and Chris Holmes.
<b>Targeting Relative Risk Heterogeneity with Causal Forests.</b> <i>In prep</i>, 2023. 
[<a href="https://drive.google.com/file/d/1lqC8FxTEPpnpy3gaVjgr9AHsfb03YWxS/view?usp=sharing">paper</a>]

Susan Athey, Julie Tibshirani and Stefan Wager.
<b>Generalized Random Forests.</b> <i>Annals of Statistics</i>, 47(2), 2019.
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>,
<a href="https://arxiv.org/abs/1610.01271">arxiv</a>]

Stefan Wager and Susan Athey.
<b>Estimation and Inference of Heterogeneous Treatment Effects using Random Forests.</b>
<i>Journal of the American Statistical Association</i>, 113(523), 2018.
[<a href="https://www.tandfonline.com/eprint/v7p66PsDhHCYiPafTJwC/full">paper</a>,
<a href="https://arxiv.org/abs/1510.04342">arxiv</a>]