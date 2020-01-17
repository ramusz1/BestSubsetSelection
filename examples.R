generate_example <- function(observations, predictors) {
  X <- matrix(rnorm(observations * predictors), byrow = TRUE, nrow = observations)
  beta <- rnorm(predictors, mean= 1)
  z <- rbinom(predictors, 1, 0.5)
  beta <- as.matrix(beta * z)
  y <- X %*% beta
  
  list(
    X = X,
    beta = beta,
    y = y,
    z = z
  )
}

generate_examples <- function(configs) {
  lapply(configs, function(config) {
    generate_example(
      observations = config$observations,
      predictors = config$predictors
    )
  })
}

examples_small_configs <- list(
  list(observations = 3, predictors = 2),
  list(observations = 5, predictors = 3),
  list(observations = 7, predictors = 5),
  list(observations = 9, predictors = 8)
)

examples_const_small_obs_configs <- list(
  list(observations = 3, predictors = 10),
  list(observations = 3, predictors = 20),
  list(observations = 3, predictors = 30),
  list(observations = 3, predictors = 40),
  list(observations = 3, predictors = 50)
)

examples_const_small_preds_configs <- list(
  list(observations = 10, predictors = 5),
  list(observations = 20, predictors = 5),
  list(observations = 30, predictors = 5),
  list(observations = 40, predictors = 5),
  list(observations = 50, predictors = 5)
)

examples_large_configs <- list(
  list(observations = 30, predictors = 30),
  list(observations = 40, predictors = 40),
  list(observations = 50, predictors = 50),
  list(observations = 60, predictors = 60)
)

examples <- list(
  small = generate_examples(examples_small_configs),
  small_obs = generate_examples(examples_const_small_obs_configs),
  small_pred = generate_examples(examples_const_small_preds_configs)
  #large= generate_examples(examples_large_configs)
)