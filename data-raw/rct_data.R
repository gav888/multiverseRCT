# Generate sample RCT data for examples
set.seed(123)
n <- 200
rct_data <- data.frame(
  id = 1:n,
  treatment_group = rep(c(0, 1), each = n/2),
  baseline = rnorm(n),
  age = sample(18:65, n, replace = TRUE)
)
rct_data$outcome_score <- 0.5 * rct_data$treatment_group +
  0.3 * rct_data$baseline + rnorm(n)

# Add binary and count outcomes for examples
rct_data$success <- as.integer(runif(n) < (0.3 + 0.2 * rct_data$treatment_group))
rct_data$count_outcome <- rpois(n, exp(1 + 0.3 * rct_data$treatment_group))

# Create data directory if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Save dataset
usethis::use_data(rct_data, overwrite = TRUE)
