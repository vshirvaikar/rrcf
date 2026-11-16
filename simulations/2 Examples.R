# Generate synthetic covariates
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

# Split training and testing data
idx <- sample.int(n, floor(0.8 * n))
X.train <- X[idx, , drop = FALSE]
Y.train <- Y[idx]
W.train <- W[idx]
X.test  <- X[-idx, , drop = FALSE]
Y.test  <- Y[-idx]
W.test  <- W[-idx]

# Next, we fit and predict all three methods on the same data. 
# This keeps comparisons fair by changing only the estimator.
# We then run a common omnibus test using the same formula.
num.trees <- 500
pval_from_rr <- function(rr.pred) {
  anova.data <- data.frame(y.test = Y.test, w.test = W.test, X.test)
  model.base <- glm(y.test ~ ., family = poisson, data = anova.data)
  model.alt <- glm(
    y.test ~ ., family = poisson,
    data = cbind(anova.data, rr.term = anova.data$w.test * log(rr.pred))
  )
  anova.out <- anova(model.base, model.alt)
  1 - pchisq(anova.out$Deviance[2], df = 1)
}

# Relative risk causal forest (RRCF)
forest.glm <- rr_causal_forest(X.train, Y.train, W.train, rct = rct_flag, num.trees = num.trees)
rr.pred.glm <- rr_predict(forest.glm, X.test)
pval.glm <- pval_from_rr(rr.pred.glm)

# Generalized random forest (GRF)
forest.grf <- causal_forest(X.train, Y.train, W.train, num.trees = num.trees)
rr.pred.grf <- rr_predict(forest.grf, X.test)
pval.grf <- pval_from_rr(rr.pred.grf)

# Model-based forest (MOB)
# MOB requires the explicit tree-control settings below 
# so its split behavior and treatment balance are stable
# and more comparable to forest-style estimation.
ctrl <- ctree_control(
  testtype = "Univ", minsplit = 2, minbucket = 2, splittry = 20,
  minprob = 0, mincriterion = 0, lookahead = TRUE, saveinfo = FALSE
)
ctrl$converged <- function(mod, data, subset) {
  all(table(data$w[subset] > 0) >= 5L)
}
train.mob <- data.frame(y = Y.train, w = W.train, X.train)
base.mob <- glm(y ~ w, data = train.mob, family = poisson())
forest.mob <- pmforest(
  base.mob, data = train.mob, ntree = num.trees,
  perturb = list(replace = FALSE, fraction = c(0.25, 0.25)),
  mtry = 2L, control = ctrl
)
coef.mob <- pmodel(x = forest.mob, newdata = X.test)
rr.pred.mob <- exp(coef.mob[, "w"])
pval.mob <- pval_from_rr(rr.pred.mob)
