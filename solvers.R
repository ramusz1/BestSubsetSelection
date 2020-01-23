cplex_solver <- function(X, beta, k) {
  run(X, beta, k)[[2]]
}

gurobi_solver <- function(X, y, k) {
  bestsubset:::bs.one.k(X, y, k, xtx = t(X) %*% X)$beta
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
  assert(dim(X)[2] >= k, "Wrong dimensions")
  object <- lars::lars(X,y,type="lasso")
  beta <- object$beta[k+1,]
  beta
}