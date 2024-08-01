#' Causal forest
#'
#' Trains a relative risk causal forest that can be used to estimate
#' relative conditional average treatment effects tau(X). When
#' the treatment assignment W is binary and unconfounded, we have
#' tau(X) = E[Y(1) | X = x] / E[Y(0) | X = x], where Y(0) and Y(1) are
#' potential outcomes corresponding to the two possible treatment states.
#'
#' @param X The covariates used in the causal regression.
#' @param Y The outcome (must be a binary numeric vector with no NAs).
#' @param W The treatment assignment (must be a binary numeric vector with no NAs).
#' @param num.trees Number of trees grown in the forest. Default is 2000.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction. Default is 0.5.
#' @param mtry Number of variables tried for each split. Default is
#'             \eqn{\sqrt p + 20} where p is the number of variables.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#'                      Default is 5.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting). Default is TRUE.
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds
#'                         to set J1 in the notation of the paper. Default is 0.5 (i.e. half of the data is used for
#'                         determining splits).
#' @param honesty.prune.leaves If TRUE, prunes the estimation sample tree such that no leaves
#'  are empty. If FALSE, keep the same tree as determined in the splits sample (if an empty leave is encountered, that
#'  tree is skipped and does not contribute to the estimate). Setting this to FALSE may improve performance on
#'  small/marginally powered data, but requires more trees. Only applies if honesty is enabled. Default is TRUE.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split. Default is 0.05.
#' @param stabilize.splits Whether or not the treatment should be taken into account when
#'                         determining the imbalance of a split. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained causal forest object. If tune.parameters is enabled,
#'  then tuning information will be included through the `tuning.output` attribute.
#'
#' @export
rr_causal_forest <- function(X, Y, W,
                             num.trees = 2000,
                             sample.fraction = 0.5,
                             mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                             min.node.size = 5,
                             honesty = TRUE,
                             honesty.fraction = 0.5,
                             honesty.prune.leaves = TRUE,
                             alpha = 0.05,
                             stabilize.splits = TRUE,
                             num.threads = NULL,
                             seed = runif(1, 0, .Machine$integer.max)) {

  has.missing.values <- validate_X(X, allow.na = TRUE)
  Y <- validate_observations(Y, X)
  W <- validate_observations(W, X)
  num.threads <- validate_num_threads(num.threads)

  model.risk = glm(Y ~ ., data = cbind(Y, X), family="poisson")
  Y.hat = predict(model.risk, X)
  Y.hat = Y.hat - min(Y.hat)
  data <- create_train_matrices(X, outcome = Y, treatment = W, sample.weights = Y.hat)

  args <- list(num.trees = num.trees,
               clusters = vector(mode = "numeric", length = 0),
               samples.per.cluster = 0,
               sample.fraction = sample.fraction,
               mtry = mtry,
               min.node.size = min.node.size,
               honesty = honesty,
               honesty.fraction = honesty.fraction,
               honesty.prune.leaves = honesty.prune.leaves,
               alpha = alpha,
               imbalance.penalty = 0,
               stabilize.splits = stabilize.splits,
               ci.group.size = 2,
               compute.oob.predictions = TRUE,
               num.threads = num.threads,
               seed = seed,
               reduced.form.weight = 0)

  forest <- do.call.rcpp(causal_train, c(data, args))
  class(forest) <- c("causal_forest", "rrcf")
  forest[["seed"]] <- seed
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["W.orig"]] <- W
  forest[["Y.hat"]] <- Y.hat
  forest[["has.missing.values"]] <- has.missing.values

  forest
}

#' Given a trained forest and test data, compute the kernel weights for each test point.
#'
#' During normal prediction, these weights (named alpha in the GRF paper) are computed as an intermediate
#' step towards producing estimates. This function allows for examining the weights directly, so they
#' could be potentially be used as the input to a different analysis.
#'
#' @param forest The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @return A sparse matrix where each row represents a test sample, and each column is a sample in the
#'         training data. The value at (i, j) gives the weight of training sample j for test sample i.
#'
#' @export
get_forest_weights <- function(forest, newdata = NULL, num.threads = NULL) {
  num.threads <- validate_num_threads(num.threads)

  forest.short <- forest[-which(names(forest) == "X.orig")]
  X <- forest[["X.orig"]]
  train.data <- create_train_matrices(X)
  args <- list(forest.object = forest.short,
               num.threads = num.threads)

  if (!is.null(newdata)) {
    test.data <- create_test_matrices(newdata)
    validate_newdata(newdata, X, allow.na = TRUE)
    do.call.rcpp(compute_weights, c(train.data, test.data, args))
  } else {
    do.call.rcpp(compute_weights_oob, c(train.data, args))
  }
}

#' Predict with a causal forest
#'
#' Gets estimates of RELATIVE tau(x) using a trained causal forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#'
#' @return Vector of RELATIVE treatment effect predictions.
#'
#' @export
rr_predict <- function(object, newdata, num.threads = NULL) {
  num.threads <- validate_num_threads(num.threads)
  X <- object[["X.orig"]]
  Y <- object[["Y.orig"]]
  W <- object[["W.orig"]]

  length = dim(newdata)[1]
  output = numeric(length)
  for(i in 1:length){
    fw = get_forest_weights(object, newdata[i,])
    df = cbind.data.frame(W, fw[1:length(Y)], fw[1:length(Y)]*Y)
    tau1 = sum(subset(df, df[,1] == 1)[,3])/sum(subset(df, df[,1] == 1)[,2])
    tau0 = sum(subset(df, df[,1] == 0)[,3])/sum(subset(df, df[,1] == 0)[,2])
    output[i] = tau1/tau0
  }
  output
}

#' ANOVA evaluation of overall detection of treatment effect heterogeneity.
#'
#' Test overall TEH detection of the forest. Conducts an ANOVA test to indicate
#' whether the predicted forest RR coefficients add any information to outcome
#' predictions, compared to baseline model with only treatment and covariates.
#'
#' @param forest The trained forest.
#' @param X The covariates used in the causal regression.
#' @param Y The outcome (must be a binary numeric vector with no NAs).
#' @param W The treatment assignment (must be a binary numeric vector with no NAs).
#'
#' @return P-value for information added by forest predictions.
#'
#' @export
rr_test_calibration <- function(forest, X, Y, W) {
  anova.data = data.frame(cbind(Y, W, X))
  predictions = rr_predict(forest, X)

  model.base = glm(Y ~ . , family = poisson, data = anova.data)
  model.forest = glm(Y ~ . , family = poisson, data = cbind(anova.data, W*log(predictions)))

  anova.test = anova(model.base, model.forest)
  result = 1-pchisq(anova.test$Deviance[2], df=1)
  result
}
