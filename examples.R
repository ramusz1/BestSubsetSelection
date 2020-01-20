generate_example <- function(observations, predictors, k) {
  X <- matrix(rnorm(observations * predictors), byrow = TRUE, nrow = observations)
  beta <- rnorm(predictors, mean= 1)
  z <- rbinom(predictors, 1, 0.5)
  beta <- as.matrix(beta * z)
  y <- X %*% beta
  
  list(
    X = X,
    beta = beta,
    k = k,
    y = y,
    z = z
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