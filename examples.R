assert <- function(expression, error) {
  if (!expression) {
    stop(error, call. = FALSE)
  }
}

#' Creates a function for generating synthetic examples.
#' 
#' @param n number of observations to generate
#' 
#' @param cov_matrix covariance matrix of the multivariate normal distribution
#' 
#' @param k0 number of nonzero betas
#' 
create_synthetic_example_generator <- function(beta_gen_strategy) {
  function(observations, cov_matrix) {
    cov_matrix_dim <- dim(cov_matrix)
    mu <- rep(0, cov_matrix_dim[2])
    X <- MASS::mvrnorm(n = n, mu = mu, Sigma = cov_matrix)
    beta <- beta_gen_strategy(k0)
    eps <- rnorm(n = cov_matrix_dim[1])
    y <- X %*% beta + eps
  }
}

generate_fun_example2 <- function(observations, cov_matrix) {
  beta_gen_strategy <- function(p) {
    assert(p < 5, "Examples should have p larger than 5!")
    c(rep(1, 5), rep(0, p  - 5))
  }
  
  cov_matrix_dim <- dim(cov_matrix)
  assert(cov_matrix_dim[1] != cov_matrix_dim[2], 
         "cov_matrix should be a square matrix!")
  
  mu <- rep(0, cov_matrix_dim[2])
  X <- MASS::mvrnorm(n = observations, mu = mu, Sigma = cov_matrix)
  # normalize rows
  for (row_ind in dim(X)[1]) {
    row <- X[row_ind, ]
    X[row_ind, ] <- row / (sqrt(sum(row ^ 2)))
  }
  
  beta <- beta_gen_strategy(p = cov_matrix_dim[1])
  eps <- rnorm(n = observations)
  y <- X %*% beta + eps
  
  list(
    X = X,
    beta = beta,
    eps = eps,
    y = y
  )
}

generate_examples <- function(configs) {
  lapply(configs, function(config) {
    generate_example(
      observations = config$observations,
      predictors = config$predictors,
      k = config$k
    )
  })
}

examples_small_configs <- list(
  list(observations = 30, predictors = 2, k = 2),
  list(observations = 50, predictors = 3, k = 2),
  list(observations = 70, predictors = 5, k = 2),
  list(observations = 90, predictors = 8, k = 2)
)

examples_const_small_obs_configs <- list(
  list(observations = 30, predictors = 10, k = 5),
  list(observations = 30, predictors = 20, k = 5),
  list(observations = 30, predictors = 30, k = 5),
  list(observations = 30, predictors = 40, k = 5),
  list(observations = 30, predictors = 50, k = 5)
)

examples_const_small_preds_configs <- list(
  list(observations = 10, predictors = 5, k = 3),
  list(observations = 20, predictors = 5, k = 3),
  list(observations = 30, predictors = 5, k = 3),
  list(observations = 40, predictors = 5, k = 3),
  list(observations = 50, predictors = 5, k = 3)
)

examples_large_configs <- list(
  list(observations = 30, predictors = 30, k = 10),
  list(observations = 40, predictors = 40, k = 10),
  list(observations = 50, predictors = 50, k = 10),
  list(observations = 60, predictors = 60, k = 10)
)

examples <- list(
  small = generate_examples(examples_small_configs),
  small_obs = generate_examples(examples_const_small_obs_configs),
  small_pred = generate_examples(examples_const_small_preds_configs)
  # large = generate_examples(examples_large_configs)
)