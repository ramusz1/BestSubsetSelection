create_solver_benchmark_fun <- function(solver_fun) {
  function(examples) {
    lapply(examples, function(example) {
      solver_fun(X = example$X, y = example$y, k = 1)
    })
  }
}

cplex_solver <- function(X, beta, k) {
  run(X, beta, k)[[2]]
}

gurobi_solver <- function(X, y, k) {
  bestsubset:::bs.one.k(X, y, k, xtx = t(X) %*% X)$beta
}

leaps_solver <- function(X, y, k) {
  regsubsets(X, y, really.big = TRUE) %>% 
    coef(id = k)
}

bs_first_order <- function(X, y, k) {
  bestsubset:::bs.proj.grad(X, y, k = k) %>% 
    as.vector()
}

#' lasso from glmnet
glmnet_solver <- function(X, y, k) {
  fit <- glmnet::glmnet(X, y)
  selections <- predict(fit, type = "nonzero")
  print(selections)
  status <- sapply(selections, function(selection){
    length(selection) == k
  })
  subset_id <- min(which(status))
  # z <- p[subset_id] use for other tests 
  coefs <- predict(fit, type = "coef")
  beta <- coefs[rownames(coefs) != '(Intercept)', subset_id]
  beta
}

lars_solver <- function(X, y, k) {
  assert(dim(X)[2] >= k, "Wrong dimensions")
  object <- lars::lars(X,y,type="lasso")
  beta <- object$beta[k+1,]
  beta
}

cplex_benchmark_fun <- create_solver_benchmark_fun(cplex_solver)
leaps_benchmark_fun <- create_solver_benchmark_fun(leaps_solver)
gurobi_benchmark_fun <- create_solver_benchmark_fun(gurobi_solver)
bs_first_order_benchmark_fun <- create_solver_benchmark_fun(bs_first_order)
lars_solver_benchmark_fun <- create_solver_benchmark_fun(lars_solver)
glmnet_solver_benchmark_fun <- create_solver_benchmark_fun(glmnet_solver)

