create_solver_benchmark_fun <- function(solver_fun) {
  function(examples) {
    lapply(examples, function(example) {
      solver_fun(X = example$X, beta = example$beta, k = example$k)
    })
  }
}

cplex_solver <- function(X, beta, k) {
  run(X, beta, k)[[2]]
}

gurobi_solver <- function(X, beta, k) {
  y <- X %*% beta
  bestsubset:::bs.one.k(X, y, k, xtx = t(X) %*% X)$beta
}

leaps_solver <- function(X, beta, k) {
  y <- X %*% beta
  regsubsets(X, y, really.big = TRUE) %>% 
    coef(id = k)
}

bs_first_order <- function(X, beta, k) {
  y <- X %*% beta
  bestsubset:::bs.proj.grad(X, y, k = k)
}

lm_solver <- function(X, beta, k) {
  y <- X %*% beta
  lm_fitted <- lm.fit(X, y)
  lm_fitted$coefficients
}


cplex_benchmark_fun <- create_solver_benchmark_fun(cplex_solver)
leaps_benchmark_fun <- create_solver_benchmark_fun(leaps_solver)
gurobi_benchmark_fun <- create_solver_benchmark_fun(gurobi_solver)
bs_first_order_benchmark_fun <- create_solver_benchmark_fun(bs_first_order)
lm_solver_benchmark_fun <- create_solver_benchmark_fun(lm_solver)
