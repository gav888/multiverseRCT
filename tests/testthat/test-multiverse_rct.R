
test_that("multiverseRCT package loads", {
  # This is a basic test just to ensure the package loads correctly
  expect_true(requireNamespace("multiverseRCT", quietly = TRUE))
})

test_that("multiverse_rct creates valid object", {
  # Create a simple example dataset
  set.seed(123)
  n <- 30
  example_data <- data.frame(
    treatment = rep(c(0, 1), each = n/2),
    covariate = rnorm(n),
    outcome = rnorm(n)
  )
  example_data$outcome <- example_data$outcome + 0.5 * example_data$treatment
  
  # Run basic multiverse analysis
  mv <- multiverse_rct(
    data = example_data,
    outcome = outcome, 
    treatment = treatment,
    preprocessing = list(none = identity),  # Only use the "none" preprocessing for testing
    covariates = ~ covariate
  )
  
  # Check that result has expected structure
  expect_s3_class(mv, "multiverse_rct")
  expect_true(nrow(mv) > 0)
  expect_true(".prep_name" %in% names(mv))
  expect_true(".model_name" %in% names(mv))
})
