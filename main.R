library(bestsubset)
library(dplyr)
library(ggplot2)
library(microbenchmark)
library(leaps)

source("./examples.R")
source("./solvers.R")

set.seed(2137)

run_benchmarks <- function(examples) {
  microbenchmark(
    cplex = cplex_benchmark_fun(examples),
    leaps = leaps_benchmark_fun(examples),
    gurobi = gurobi_benchmark_fun(examples),
    bs_first_order = bs_first_order_benchmark_fun(examples),
    lm = lm_solver_benchmark_fun(examples),
    times = 10
  )
}

mse_error <- function(example, result) {
  y_result <- example$X %*% result
  mean((example$y - y_result)^2)
}

mse_solvers <- list(
  cplex = cplex_benchmark_fun,
  gurobi = gurobi_benchmark_fun,
  bs_first_order = bs_first_order_benchmark_fun,
  lm = lm_solver_benchmark_fun
)

xor_error <- function(example, result) {
  z <- result != 0
  sum(xor(example$z, result))
}

xor_solvers <- list(
  cplex = cplex_benchmark_fun,
  gurobi = gurobi_benchmark_fun,
  bs_first_order = bs_first_order_benchmark_fun,
  leaps = leaps_benchmark_fun
)

run_precision_benchmarks <- function(examples, solvers, error_fun) {
  run_single_precision_benchmarks <- function(solver_benchmark_fun) {
    results <- solver_benchmark_fun(examples)
    errors <- sapply(seq_along(results), function(k){
      error_fun(examples[[k]], results[[k]])
    })
    mean(errors)
  }
  
  sapply(solvers, run_single_precision_benchmarks)
}

visualize_benchmarks <- function(benchmark) {
  df <- data.frame(
    expr = benchmark$expr,
    time = benchmark$time,
    stringsAsFactors = FALSE
  )
  
  ordered_df <- df %>% 
    group_by(expr) %>% 
    summarise(mean_time = mean(time)) %>% 
    arrange(desc(mean_time))
  
  expr_levels <- ordered_df$expr
  ordered_df$expr <- factor(ordered_df$expr, rev(expr_levels))
  ordered_df %>% 
    ggplot(aes(x = expr, y = mean_time, fill = expr)) + 
    geom_col() + 
    labs(x = "Solver", y = "Mean time", fill="Solver") + 
    coord_flip()
}

visualize_error <- function(error_benchmarks, error_name) {
  solvers <- names(error_benchmarks)
  names(error_benchmarks) <- NULL
  
  df <- data.frame(
    solver = solvers,
    error = error_benchmarks,
    stringsAsFactors = FALSE
  )
  
  solver_levels <- df %>% 
    arrange(error) %>% 
    pull(solver)
  
  df$solver <- factor(df$solver, solver_levels)
  
  df %>% 
    ggplot(aes(x = solver, y = error, fill = solver)) + 
    geom_col() + 
    labs(x = "Solver", y = error_name, fill="Solver") + 
    coord_flip()
}

# Run benchmarks
benchmarks <- lapply(examples, run_benchmarks)
mse_benchmarks <- lapply(examples, function(example) run_precision_benchmarks(example, mse_solvers, mse_error))
xor_benchmarks <- lapply(examples, function(example) run_precision_benchmarks(example, xor_solvers, xor_error))

# Visualize benchmarks
visualizations <- lapply(benchmarks, visualize_benchmarks)
mse_visualizations <- lapply(mse_benchmarks, function(mse_benchmark) {
  visualize_error(mse_benchmark, error_name = "mse error")
})

xor_visualizations <- lapply(xor_benchmarks, function(benchmark) {
  visualize_error(benchmark, error_name = "xor error")
})
