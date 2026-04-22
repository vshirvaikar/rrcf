library(model4you)

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
X.test <- X[-idx, , drop = FALSE]
Y.test <- Y[-idx]
W.test <- W[-idx]

# Fit relative risk model-based forest
mob.train <- data.frame(y = Y.train, w = W.train, X.train)
mob.base <- glm(y ~ w, data = mob.train, family = binomial)
forest.rr <- pmforest(mob.base, data = mob.train)

# Predict relative conditional treatment effect
mob.test.1 <- data.frame(w = 1, X.test)
mob.test.0 <- data.frame(w = 0, X.test)
pred.1 <- predict(forest.rr, mob.test.1)[[1]]
pred.0 <- predict(forest.rr, mob.test.0)[[1]]
rr.pred <- unname(pred.1 / pred.0)

# Test for treatment effect heterogeneity
anova.data <- data.frame(y.test = Y.test, w.test = W.test, X.test)
model.base <- glm(y.test ~ ., family = poisson, data = anova.data)
model.mob <- glm(
  y.test ~ ., family = poisson,
  data = cbind(anova.data, rr.mob = anova.data$w.test * log(rr.pred))
)
anova.mob <- anova(model.base, model.mob)
pval <- 1 - pchisq(anova.mob$Deviance[2], df = 1)