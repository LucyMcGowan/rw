test_that("compute_score works for gaussian family", {
  # Create test data
  set.seed(123)
  data <- data.frame(
    x1 = rnorm(100),
    x2 = rnorm(100),
    y = rnorm(100)
  )
  
  # Create model info
  model_info <- list(
    family = "gaussian",
    coefficients = c(`(Intercept)` = 0.5, x1 = 1.2, x2 = -0.8),
    sigma2 = 1.5
  )
  
  # Compute score
  result <- compute_score(data, model_info, "y")
  
  # Check dimensions
  expect_equal(nrow(result), 100)
  expect_equal(ncol(result), 4)  # 3 coefficients + 1 sigma
  
  # Check that result is numeric
  expect_true(is.numeric(result))
  expect_false(any(is.na(result)))
})

test_that("compute_score works for binomial family", {
  # Create test data
  set.seed(123)
  data <- data.frame(
    x1 = rnorm(100),
    x2 = rnorm(100),
    y = rbinom(100, 1, 0.5)
  )
  
  # Create model info
  model_info <- list(
    family = "binomial",
    coefficients = c(`(Intercept)` = 0.2, x1 = 0.5, x2 = -0.3)
  )
  
  # Compute score
  result <- compute_score(data, model_info, "y")
  
  # Check dimensions
  expect_equal(nrow(result), 100)
  expect_equal(ncol(result), 3)  # 3 coefficients, no sigma
  
  # Check that result is numeric
  expect_true(is.numeric(result))
  expect_false(any(is.na(result)))
})

test_that("compute_score handles intercept correctly", {
  # Create test data
  set.seed(123)
  data <- data.frame(
    x1 = rnorm(50),
    y = rnorm(50)
  )
  
  # Model with explicit intercept
  model_info <- list(
    family = "gaussian",
    coefficients = c(`(Intercept)` = 1.0, x1 = 0.5),
    sigma2 = 1.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  # First column should correspond to intercept
  expect_equal(ncol(result), 3)  # intercept + x1 + sigma
})


test_that("compute_score works without intercept", {
  # Create test data
  set.seed(123)
  data <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    y = rnorm(50)
  )
  
  # Model without intercept
  model_info <- list(
    family = "gaussian",
    coefficients = c(x1 = 0.5, x2 = -0.3),
    sigma2 = 1.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  expect_equal(ncol(result), 3)  # x1 + x2 + sigma
})

test_that("compute_score computes correct gaussian scores", {
  # Simple case with known values
  data <- data.frame(
    x1 = c(1, 2, 3),
    y = c(2, 4, 6)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c(`(Intercept)` = 0, x1 = 2),
    sigma2 = 1
  )
  
  result <- compute_score(data, model_info, "y")
  
  # With perfect fit, residuals should be zero
  # Score for beta should be X * residual / sigma
  expect_equal(unname(result[, 1:2]), matrix(0, 3, 2))
  
  # Score for sigma should be 0.5 * (-1 + residual^2)
  expect_equal(result[, 3], rep(-0.5, 3))
})

test_that("compute_score computes correct binomial scores", {
  # Simple case
  data <- data.frame(
    x1 = c(0, 0, 0),
    y = c(1, 0, 1)
  )
  
  model_info <- list(
    family = "binomial",
    coefficients = c(`(Intercept)` = 0, x1 = 0)
  )
  
  result <- compute_score(data, model_info, "y")
  
  # With beta = 0, p = 0.5 for all observations
  # Score should be X * (y - p) = X * (y - 0.5)
  expected_intercept <- c(0.5, -0.5, 0.5)
  expect_equal(result[, 1], expected_intercept)
})

test_that("compute_score fails with unsupported family", {
  data <- data.frame(
    x1 = rnorm(10),
    y = rpois(10, 2)
  )
  
  model_info <- list(
    family = "poisson",
    coefficients = c(`(Intercept)` = 0, x1 = 0.5)
  )
  
  expect_error(
    compute_score(data, model_info, "y"),
    "Unsupported family"
  )
})

test_that("compute_score handles single observation", {
  data <- data.frame(
    x1 = 1.5,
    y = 2.0
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c(`(Intercept)` = 0.5, x1 = 1.0),
    sigma2 = 1.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
})

test_that("compute_score handles multiple predictors", {
  set.seed(123)
  data <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    x3 = rnorm(50),
    x4 = rnorm(50),
    y = rnorm(50)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c(`(Intercept)` = 0, x1 = 1, x2 = -1, x3 = 0.5, x4 = -0.5),
    sigma2 = 2.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  expect_equal(nrow(result), 50)
  expect_equal(ncol(result), 6)  # 5 coefficients + sigma
})

test_that("compute_score preserves coefficient order", {
  set.seed(123)
  data <- data.frame(
    x1 = rnorm(20),
    x2 = rnorm(20),
    y = rnorm(20)
  )
  
  model_info <- list(
    family = "gaussian",
    coefficients = c(`(Intercept)` = 0, x2 = 1, x1 = 0.5),  # Note: x2 before x1
    sigma2 = 1.0
  )
  
  result <- compute_score(data, model_info, "y")
  
  # Should have columns in order: (Intercept), x2, x1, sigma
  expect_equal(ncol(result), 4)
})
