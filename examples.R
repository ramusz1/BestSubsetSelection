library(magrittr)

assert <- function(expression, error) {
  if (!expression) {
    stop(error, call. = FALSE)
  }
}


normalize_x <- function(X) {
  m <- dim(X)[2]
  for (col_ind in seq_len(m)) {
    col <- X[, col_ind]
    X[, col_ind] <- col / (sqrt(sum(col ^ 2)))
  }
  X
}

get_diabetes_data <- function() {
  data(diabetes, package = "lars")
  X <- diabetes$x2
  # unit l2
  for (col_ind in seq_len(dim(X)[2])) {
    # standardize to zero mean
    col <- X[, col_ind] - mean(X[, col_ind])
    X[, col_ind] <- col / (sqrt(sum(col ^ 2)))
  }
  y <- diabetes$y
  
  list(
    X = X,
    y = y
  )
}

#' Noise genarator based on signal noise ratio
get_noise_from_snr <- function(y, snr){
  sd <- sqrt(var(y) / snr)
  rnorm(n = length(y), sd=sd)
}

create_identity_matrix <- function(p) {
  diag(p)
}

#' Creates a function for generating synthetic examples.
#' 
#' @param k0 number of nonzero betas
#' 
create_synthetic_example_generator <- function(beta_gen_strategy, cov_matrix_gen_strategy) {
  function(observations, p, snr=7) {
    cov_matrix <- cov_matrix_gen_strategy(p)
    cov_matrix_dim <- dim(cov_matrix)
    assert(cov_matrix_dim[1] == cov_matrix_dim[2], 
           "cov_matrix should be a square matrix!")
    
    mu <- rep(0, cov_matrix_dim[2])
    X <- MASS::mvrnorm(n = observations, mu = mu, Sigma = cov_matrix)
    # normalize cols
    X <- normalize_x(X)
    
    beta <- beta_gen_strategy(p = cov_matrix_dim[1])
    noiseless_y <- X %*% beta
    eps <- get_noise_from_snr(noiseless_y, snr)
    y <- noiseless_y + eps
    
    list(
      X = X,
      beta = beta,
      eps = eps,
      y = y
    )
  }
}



#' Example 1 generation
create_cov_matrix_gen_strategy <- function(ro) {
  function(p) {
    cov_matrix <- matrix(rep(0, p * p), nrow = p, ncol = p)
    for (row in seq_len(p)) {
      for (col in seq_len (p)) {
        cov_matrix[row, col] <- abs(row - col) 
      }
    }
    
    ro ^ cov_matrix
  }
}

beta_gen_strategy_example1 <- function(p) {
  assert(p > 10, "Examples should have p larger than 10!")
  k0 <- 10
  c(rep(0, p-k0), rep(1, k0))
}

#' Example 2 generation
beta_gen_strategy_example2 <- function(p) {
  assert(p > 5, "Examples should have p larger than 5!")
  c(rep(1, 5), rep(0, p  - 5))
}
generate_fun_example2 <- create_synthetic_example_generator(
  beta_gen_strategy = beta_gen_strategy_example2,
  cov_matrix_gen_strategy = create_identity_matrix
)

#' Example 3 generation
beta_gen_strategy_example3 <- function(p) {
  assert(p > 10, "Examples should have p larger than 10!")
  c(0.5 + 9.5 * ((seq_len(10) - 1) / 10), rep(0, p - 10))
}

#' Example 4 generation
beta_gen_strategy_example4 <- function(p) {
  assert(p > 6, "Examples should have p larger than 10!")
  c(-10, -6, -2, 2, 6, 10, rep(0, p - 6))
}
generate_fun_example4 <- create_synthetic_example_generator(
  beta_gen_strategy = beta_gen_strategy_example4,
  cov_matrix_gen_strategy = create_identity_matrix
)


## Examples generaton
create_problem_examples <- function(configs, example_generator) {
  generate_problem_example <- function(config) {
    example <- example_generator(
      observations = config$observations,
      p = config$p
    )
    list(
      problem = example,
      k = config$k
    )
  }
  lapply(configs, generate_problem_example)
}

create_problem_examples_like_example1 <- function(configs) {
  generate_problem_example <- function(config) {
    cov_matrix_gen_strategy <- create_cov_matrix_gen_strategy(ro = config$ro)
    example_generator <- create_synthetic_example_generator(
      beta_gen_strategy = beta_gen_strategy_example1,
      cov_matrix_gen_strategy = cov_matrix_gen_strategy
    )
    example <- example_generator(
      observations = config$observations,
      p = config$p,
      snr = config$snr
    )
    list(
      problem = example,
      k = config$k,
      snr = config$snr,
      ro = config$ro
    )
  }
  lapply(configs, generate_problem_example)
}

fixed_variables_changing_observations_example_config <- list(
  "500 observations" = list(observations = 500, p = 10, k = 7),
  "1500 observations" = list(observations = 1500, p = 10, k = 7),
  "2500 observations" = list(observations = 2500, p = 10, k = 7),
  "3500 observations" = list(observations = 3500, p = 10, k = 7),
  "4500 observations" = list(observations = 4500, p = 10, k = 7)
)

fixed_observations_changing_variables_example_config <- list(
  "10 variables" = list(observations = 5000, p = 10, k = 5),
  "20 variables" = list(observations = 5000, p = 20, k = 10),
  "30 variables" = list(observations = 5000, p = 30, k = 15)
)

precision_and_best_subset_examples_config <- list(
  "ro=0.5, snr=1.58" = list(observations = 500, p = 100, k = 100, ro = 0.5, snr = 1.58),
  "ro=0.5, snr=3.17" = list(observations = 500, p = 100, k = 100, ro = 0.5, snr = 3.17),
  "ro=0.5, snr=6.33" = list(observations = 500, p = 100, k = 100, ro = 0.5, snr = 6.33),
  "ro=0.8, snr=1.74" = list(observations = 500, p = 100, k = 100, ro = 0.8, snr = 1.74),
  "ro=0.8, snr=3.48" = list(observations = 500, p = 100, k = 100, ro = 0.8, snr = 3.48),
  "ro=0.8, snr=6.97" = list(observations = 500, p = 100, k = 100, ro = 0.8, snr = 6.97),
  "ro=0.9, snr=2.18" = list(observations = 500, p = 100, k = 100, ro = 0.9, snr = 2.18),
  "ro=0.9, snr=4.37" = list(observations = 500, p = 100, k = 100, ro = 0.9, snr = 4.37),
  "ro=0.9, snr=8.73" = list(observations = 500, p = 100, k = 100, ro = 0.9, snr = 8.73)
)

examples <- list(
  fixed_variables_changing_observations = create_problem_examples(
    configs = fixed_variables_changing_observations_example_config,
    example_generator = generate_fun_example2
  ),
  fixed_observations_changing_variables = create_problem_examples(
    configs = fixed_observations_changing_variables_example_config,
    example_generator = generate_fun_example2
  ),
  precision_and_best_subset_exmaple = create_problem_examples_like_example1(
    configs = precision_and_best_subset_examples_config
  )
)
