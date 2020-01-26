
#' Create CPLEX MIQP solver
#' @param start can be one of: 
#'        'warm' - starting using first order solution
#'        'cold' - solver starts without any priors 
#'        'mild' - starting point based on lm 
#'        'theory' - using theoretical bounds and random starting point
get_cplex_solver <- function(start) {
  function(X, y, k) {
    bestsubset:::run_bs(X, y, k, start, run_cplex_miqp)
  }
}

#' Create GUROBI MIQP solver
#' @param start can be one of: 
#'        'warm' - starting using first order solution
#'        'cold' - solver starts without any priors 
#'        'mild' - starting point based on lm 
#'        'theory' - using theoretical bounds and random starting point
get_gurobi_solver <- function(start) {
  function(X, y, k){
    bestsubset:::run_bs(X, y, k, start, bestsubset:::miqp_bs)$beta
  }
}

leaps_solver <- function(X, y, k) {
  leaps::leaps(X, y) %>% 
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
  coefs <- predict(fit, type = "coef")
  beta <- coefs[rownames(coefs) != '(Intercept)', subset_id]
  beta
}

lars_solver <- function(X, y, k) {
  object <- lars::lars(X,y,type="lasso")
  beta <- object$beta[k+1,]
  beta
}