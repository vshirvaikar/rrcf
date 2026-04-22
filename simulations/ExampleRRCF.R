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
