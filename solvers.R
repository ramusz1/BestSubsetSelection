create_solver_benchmark_fun <- function(solver_fun) {
  function(examples) {
    lapply(examples, function(example) {
      solver_fun(X = example$X, beta = example$beta)
    })
  }
}

cplex_solver <- function(X, beta) {
  run(X, beta)[[2]]
}

gurobi_solver <- function(X, beta) {
  y <- X %*% beta
  bs(X, y)
}

leaps_solver <- function(X, beta) {
  y <- X %*% beta
  regsubsets(X, y, really.big = TRUE)
}


cplex_benchmark_fun <- create_solver_benchmark_fun(cplex_solver)
leaps_benchmark_fun <- create_solver_benchmark_fun(leaps_solver)
gurobi_benchmark_fun <- create_solver_benchmark_fun(gurobi_solver)