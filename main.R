library(bestsubset)
library(dplyr)
library(ggplot2)
library(microbenchmark)
library(leaps)
library(parallel)

set.seed(2137)

source("./examples.R")
source("./performance.R")
source("./solvers.R")
source("./visualizations.R")



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
  
  names(exprs) <- names(solvers)
  microbenchmark(
    list = exprs,
    times = times
  )
}


solvers_performance_benchmark <- function(solvers, problem_examples, cl) {
  cur_env <- environment()
  clusterExport(cl = cl, "solvers", envir = cur_env)
  clusterExport(cl = cl, "problem_examples", envir = cur_env)
  clusterExport(cl = cl, "calculate_prediction_performance", envir = cur_env)
  Sys.sleep(10)
  
  single_solver_performance_benchmark <- function(solvers, problem_example) {
    single_solver_benchmark <- function(solver) {
      print("Solver!")
      k_list <- seq(1, problem_example$k)
      model_results <- parLapply(cl = cl, k_list, function(k) {
        library(bestsubset)
        library(magrittr)
        result <- solver(X = problem_example$problem$X,
                         y = problem_example$problem$y,
                         k = k)
        performance <- calculate_prediction_performance(
          X = problem_example$problem$X,
          beta = problem_example$problem$beta,
          beta_hat = result
        )
        
        list(
          k = k,
          performance = performance
        )
      }) %>% bind_rows()
      
      best_result <- model_results %>% 
        filter(performance == min(performance)) %>% 
        filter(k == min(k))
      
      data.frame(
        nonzeros = best_result$k,
        performance = best_result$performance,
        snr = problem_example$snr,
        ro = problem_example$ro,
        stringsAsFactors = FALSE
      )
    }
    df <- lapply(solvers, single_solver_benchmark) %>% 
      bind_rows()
    df$solver <- names(solvers)
    df
  }
  
  bind_rows(lapply(problem_examples, function(problem_example) {
    single_solver_performance_benchmark(solvers, problem_example)
  }))
}


# Time benchmark for growing number of observations
solvers <- list(
  cplex_warm = get_cplex_solver('warm'),
  cplex_cold = get_cplex_solver('cold'),
  gurobi_warm = get_gurobi_solver('warm'),
  gurobi_cold = get_gurobi_solver('cold'),
  first_order = bs_first_order,
  leaps = leaps_solver,
  lars = lars_solver
)
fixed_variables_changing_observations_benchmarks <- lapply(solvers, function(solver) {
  solver_time_benchmark(solver, examples$fixed_variables_changing_observations, 100)
})
fixed_variables_changing_observations_plot <- plot_solver_times(
  benchmarks = fixed_variables_changing_observations_benchmarks
)

# Time benchmark for growing number of variables
solvers <- list(
  cplex_warm = get_cplex_solver('warm'),
  cplex_cold = get_cplex_solver('cold'),
  gurobi_warm = get_gurobi_solver('warm'),
  gurobi_cold = get_gurobi_solver('cold'),
  first_order = bs_first_order,
  lars = lars_solver
)


fixed_observations_changing_variables_benchmarks <- lapply(solvers, function(solver) {
  solver_time_benchmark(solver, examples$fixed_observations_changing_variables, 100)
})
fixed_observations_changing_variables_plot <- plot_solver_times(
  benchmarks = fixed_observations_changing_variables_benchmarks
)

# Precision and subset size benchmark
num_cores <- detectCores()
cl <- makeCluster(num_cores)
solvers <- list(
  cplex_warm = get_cplex_solver('warm'),
  cplex_cold = get_cplex_solver('cold'),
  gurobi_warm = get_gurobi_solver('warm'),
  gurobi_cold = get_gurobi_solver('cold'),
  first_order = bs_first_order,
  lars = lars_solver
)

precision_and_best_subset_benchmarks <- solvers_performance_benchmark(
  solvers = solvers,
  problem_examples = examples$precision_and_best_subset_exmaple,
  cl = cl
)
#saveRDS(precision_and_best_subset_benchmarks, "backup-full.RDS")
plot_performance_benchmark(precision_and_best_subset_benchmarks)
cl <- NULL















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

# TODO update
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