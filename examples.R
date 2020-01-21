assert <- function(expression, error) {
  if (!expression) {
    stop(error, call. = FALSE)
  }
}


#' Noise genarator based on signal noise ratio
get_noise_from_snr <- function(y, snr){
  sd <- sqrt(var(y) / snr)
  rnorm(n = length(y), sd=sd)
}

# Check if snr noise generator works
# X <- 1:100
# b <- 1
# noiseless_y <- X * b
# snrs <- c(2, 4, 8)
# y_with_noise <- lapply(snrs, function(snr){
#   noise <- get_noise_from_snr(noiseless_y, snr)
#   data.frame(y=noiseless_y + noise, snr=snr, x=X)
# })
# df <- do.call(rbind, y_with_noise)
# df <- rbind(df, data.frame(y=noiseless_y, snr='inf', x=X))
# library(ggplot2)
# ggplot(df, aes(x=x, y=y, color=snr)) +
#   geom_point()

#' Creates a function for generating synthetic examples.
#' 
#' @param n number of observations to generate
#' 
#' @param cov_matrix covariance matrix of the multivariate normal distribution
#' 
#' @param k0 number of nonzero betas
#' 
create_synthetic_example_generator <- function(beta_gen_strategy) {
  function(observations, cov_matrix, snr=7) {
    cov_matrix_dim <- dim(cov_matrix)
    assert(cov_matrix_dim[1] == cov_matrix_dim[2], 
           "cov_matrix should be a square matrix!")
    
    mu <- rep(0, cov_matrix_dim[2])
    X <- MASS::mvrnorm(n = observations, mu = mu, Sigma = cov_matrix)
    # normalize rows
    for (row_ind in seq_len(dim(X)[1])) {
      row <- X[row_ind, ]
      X[row_ind, ] <- row / (sqrt(sum(row ^ 2)))
    }
    
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
create_cov_matrix_for_example1 <- function(ro, p) {
  cov_matrix <- matrix(rep(0, p * p), nrow = p, ncol = p)
  for (row in seq_len(p)) {
    for (col in seq_len (p)) {
      cov_matrix[row, col] <- abs(row - col) 
    }
  }
  
  ro ^ cov_matrix
}

beta_gen_strategy_example1 <- function(p) {
  assert(p > 10, "Examples should have p larger than 10!")
  k0 <- round(runif(1, min = 5, max = 10))
  c(rep(0, k0), rep(1, p - k0))
}
generate_fun_example1 <- create_synthetic_example_generator(
  beta_gen_strategy = beta_gen_strategy_example1)
example1_configs <- list(
  list(observations = 10, ro = 0.7, p = 11, snr = 3),
  list(observations = 10, ro = 0.7, p = 11, snr = 7),
  list(observations = 10, ro = 0.7, p = 12, snr = 3),
  list(observations = 10, ro = 0.7, p = 12, snr = 7),
  list(observations = 10, ro = 0.7, p = 13, snr = 7)
)

example1_cases <- lapply(example1_configs, function(config) {
  cov_matrix <- create_cov_matrix_for_example1(ro = config$ro, p = config$p)
  generate_fun_example1(
    observations = config$observations,
    cov_matrix = cov_matrix,
    snr=config$snr
  )
})

#' Example 2 generation
beta_gen_strategy_example2 <- function(p) {
  assert(p > 5, "Examples should have p larger than 5!")
  c(rep(1, 5), rep(0, p  - 5))
}
generate_fun_example2 <- create_synthetic_example_generator(beta_gen_strategy = beta_gen_strategy_example2)
example2_configs <- list(
  list(observations = 10, p = 11),
  list(observations = 10, p = 12),
  list(observations = 10, p = 13),
  list(observations = 10, p = 10),
  list(observations = 10, p = 10)
)
example2_cases <- lapply(example1_configs, function(config) {
  cov_matrix <- diag(config$p)
  generate_fun_example2(
    observations = config$observations,
    cov_matrix = cov_matrix
  )
})


#' Example 3 generation
beta_gen_strategy_example3 <- function(p) {
  assert(p > 10, "Examples should have p larger than 10!")
  c(0.5 + 9.5 * ((seq_len(10) - 1) / 10), rep(0, p - 10))
}
generate_fun_example3 <- create_synthetic_example_generator(beta_gen_strategy = beta_gen_strategy_example3)
example3_configs <- list(
  list(observations = 10, p = 11),
  list(observations = 10, p = 12),
  list(observations = 10, p = 13),
  list(observations = 10, p = 10),
  list(observations = 10, p = 10)
)
example3_cases <- lapply(example1_configs, function(config) {
  cov_matrix <- diag(config$p)
  generate_fun_example3(
    observations = config$observations,
    cov_matrix = cov_matrix
  )
})


#' Example 4 generation
beta_gen_strategy_example4 <- function(p) {
  assert(p > 6, "Examples should have p larger than 10!")
  c(-10, -6, -2, 2, 6, 10, rep(0, p - 6))
}
generate_fun_example4 <- create_synthetic_example_generator(beta_gen_strategy = beta_gen_strategy_example4)
example4_configs <- list(
  list(observations = 10, p = 11),
  list(observations = 10, p = 12),
  list(observations = 10, p = 13)
)
example4_cases <- lapply(example1_configs, function(config) {
  cov_matrix <- diag(config$p)
  generate_fun_example4(
    observations = config$observations,
    cov_matrix = cov_matrix
  )
})  

examples <- list(
  example1 = example1_cases,
  example2 = example2_cases,
  example3 = example3_cases,
  example4 = example4_cases
)