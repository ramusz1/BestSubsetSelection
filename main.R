library(bestsubset)
library(dplyr)
library(ggplot2)
library(microbenchmark)
library(leaps)

set.seed(2137)

source("./examples.R")
source("./solvers.R")

solver_time_benchmark <- function(solver, problem_examples, times) {
  exprs <- lapply(problem_examples, function(problem_example) {
    rlang::expr(solver(X = !!(problem_example$problem$X), y = !!(problem_example$problem$y), k = !!(problem_example$k)))
  })
  
  microbenchmark(
    list = exprs,
    times = times
  )
}

solver_comparison_time_benchmark <- function(solvers, problem_example, times) {
  exprs <- lapply(names(solvers), function(solver_name) {
    rlang::expr(solvers[[!!solver_name]](X = !!(problem_example$problem$X), y = !!(problem_example$problem$y), k = !!(problem_example$k)))
  })
  
  microbenchmark(
    list = exprs,
    times = times
  )
}

# Single solver many problems example
k_diff_problem_examples <- lapply(3:8, function(k) {
  list(
    problem = example1_cases[[1]],
    k = k
  )
})
names(k_diff_problem_examples) <- paste("Case: k =", 3:8) # Must be named!
solver_time_benchmark(gurobi_solver, k_diff_problem_examples, 4)

# Many solvers single problem example
solvers <- list(
  leaps = leaps_solver,
  gurobi = gurobi_solver
)
problem_example <- list(
  problem = example1_cases[[1]],
  k = 10
)

sample_comparison <- solver_comparison_time_benchmark(solvers = solvers, problem_example, 5)
boxplot(sample_comparison)

# Different Cplex configs on single problem example 
solvers <- list(
  cplex_warm = get_cplex_solver('warm'),
  cplex_mild = get_cplex_solver('mild'),
  cplex_cold = get_cplex_solver('cold'),
  cplex_theory = get_cplex_solver('theory') # fails
)

problem_example <- list(
  problem = example5_cases[[1]],
  k = 5
)

reticulate::source_python("main.py")
sample_comparison <- solver_comparison_time_benchmark(solvers = solvers, problem_example, 3)
boxplot(sample_comparison)


mse_error <- function(example, result) {
  y_result <- example$X %*% result
  mean((example$y - y_result)^2)
}

mse_solvers <- list(
  gurobi = gurobi_benchmark_fun,
  bs_first_order = bs_first_order_benchmark_fun
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

# Save plots
fig_dir <- "./figures/"
lapply(names(visualizations), function(n) {
  filename <- paste0("time_benchmark-", n, ".pdf")
  ggsave(file.path(fig_dir, filename), visualizations[[n]])
})

lapply(names(mse_visualizations), function(n) {
  filename <- paste0("mse_benchmark-", n, ".pdf")
  ggsave(file.path(fig_dir, filename), mse_visualizations[[n]])
})

lapply(names(xor_visualizations), function(n) {
  filename <- paste0("xor_benchmark-", n, ".pdf")
  ggsave(file.path(fig_dir, filename), xor_visualizations[[n]])
})