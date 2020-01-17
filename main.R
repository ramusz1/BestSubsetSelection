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
    times = 10
  )
}

run_precision_benchmarks <- function(examples) {
  run_single_precision_benchmarks <- function(solver_benchmark_fun) {
    calculate_error <- function(example, result) {
      y_result <- example$X %*% result
      (example$y - y_result)^2
    }
    
    results <- solver_benchmark_fun(examples)
    errors <- sapply(seq_along(results), function(k){
      calculate_error(examples[[k]], results[[k]])
    })
    
    mean(errors)
  }
  
  solvers <- list(
    cplex = cplex_benchmark_fun,
    leaps = leaps_benchmark_fun,
    gurobi = gurobi_benchmark_fun
  )
  
  lapply(solvers, run_single_precision_benchmarks)
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
    labs(y = "Mean time", fill="Solver") + 
    coord_flip()
}

# Run benchmarks
benchmarks <- lapply(examples, run_benchmarks)
precision_benchmarks <- lapply(examples, run_precision_benchmarks)

# Visualize benchmarks
visualizations <- lapply(benchmarks, visualize_benchmarks)


