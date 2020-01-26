calculate_prediction_performance <- function(X, beta, beta_hat) {
  l2_squared <- function(x) {
    sum(x ^ 2)
  }
  l2_squared(X %*% beta - X%*% beta_hat) / l2_squared(X %*% beta)
}
