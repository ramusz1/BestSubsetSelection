library(dplyr)
library(ggplot2)
library(patchwork)

plot_solver_times <- function(benchmarks) {
  benchmark_dfs <- lapply(names(benchmarks), function(solver_name) {
    benchmark <- benchmarks[[solver_name]]
    df <- data.frame(
      expr = benchmark$expr,
      time = benchmark$time,
      stringsAsFactors = FALSE
    )
    df$solver <- solver_name
    df %>% group_by(expr, solver) %>% summarize(mean_time = mean(time))
  })
  df <- bind_rows(benchmark_dfs)
  ggplot(df, aes(x = expr, y = solver, fill = mean_time)) + geom_tile() + scale_fill_viridis_c()
}

plot_performance_benchmark <- function(benchmark) {
  precision_plot <- benchmark %>% 
    mutate(solver = factor(solver),
           snr = factor(snr)) %>% 
    ggplot(aes(x = snr, y = performance)) +
    geom_bar(aes(fill = solver), position = "dodge", stat = "identity") +
    facet_wrap(. ~ ro, scales = "free")
  
  nonzeros_plot <- benchmark %>% 
    mutate(solver = factor(solver),
           snr = factor(snr)) %>% 
    ggplot(aes(x = snr, y = nonzeros)) +
    geom_bar(aes(fill = solver), position = "dodge", stat = "identity") +
    facet_wrap(. ~ ro, scales = "free")
  
  nonzeros_plot / precision_plot 
}
