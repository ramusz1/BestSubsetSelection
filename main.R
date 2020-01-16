library(reticulate)
use_python("venv/bin/python")

source_python("main.py")

sample_matrix <- matrix(
  data = runif(10 * 10),
  nrow = 10,
  byrow = TRUE
)

sample_beta <- matrix(
  data = runif(10),
  nrow = 10,
  byrow = TRUE
)

sample_result <- run(sample_matrix, sample_beta)

