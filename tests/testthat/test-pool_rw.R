test_that("pool_rw requires with_rw object", {
  expect_error(
    pool_rw(list(results = NULL))
  )
  
  expect_error(
    pool_rw(data.frame(x = 1:5))
  )
})

test_that("pool_rw returns correct structure", {
  # Mock a simple with_rw object
  mock_model1 <- list(coefficients = c(`(Intercept)` = 1.0, x1 = 0.5))
  class(mock_model1) <- "lm"
  
  mock_model2 <- list(coefficients = c(`(Intercept)` = 1.1, x1 = 0.6))
  class(mock_model2) <- "lm"
  
  with_rw_obj <- structure(
    list(
      m = 2,
      n = 100,
      results = list(
        list(
          model = mock_model1,
          U = matrix(rnorm(200), 100, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(rnorm(100), 100, 1),
          d = matrix(rnorm(100), 100, 1)
        ),
        list(
          model = mock_model2,
          U = matrix(rnorm(200), 100, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(rnorm(100), 100, 1),
          d = matrix(rnorm(100), 100, 1)
        )
      )
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  expect_s3_class(result, "pool_rw")
  expect_named(result, c("pooled", "variance", "m", "n", "call"))
  expect_equal(result$m, 2)
  expect_equal(result$n, 100)
})

test_that("pool_rw pooled estimates are average of coefficients", {
  # Create mock with_rw object
  mock_model1 <- list(coefficients = c(`(Intercept)` = 1.0, x1 = 2.0))
  class(mock_model1) <- "lm"
  
  mock_model2 <- list(coefficients = c(`(Intercept)` = 1.2, x1 = 2.4))
  class(mock_model2) <- "lm"
  
  mock_model3 <- list(coefficients = c(`(Intercept)` = 0.8, x1 = 1.6))
  class(mock_model3) <- "lm"
  
  with_rw_obj <- structure(
    list(
      m = 3,
      n = 100,
      results = list(
        list(
          model = mock_model1,
          U = matrix(0, 100, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(0, 100, 1),
          d = matrix(0, 100, 1)
        ),
        list(
          model = mock_model2,
          U = matrix(0, 100, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(0, 100, 1),
          d = matrix(0, 100, 1)
        ),
        list(
          model = mock_model3,
          U = matrix(0, 100, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(0, 100, 1),
          d = matrix(0, 100, 1)
        )
      )
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  # Expected means
  expect_equal(result$pooled$estimate[1], 1.0)  # (1.0 + 1.2 + 0.8) / 3
  expect_equal(result$pooled$estimate[2], 2.0)  # (2.0 + 2.4 + 1.6) / 3
})

test_that("pool_rw computes variance matrix correctly", {
  set.seed(123)
  
  mock_model <- list(coefficients = c(`(Intercept)` = 1.0, x1 = 0.5))
  class(mock_model) <- "lm"
  
  n <- 100
  with_rw_obj <- structure(
    list(
      m = 5,
      n = n,
      results = lapply(1:5, function(i) {
        list(
          model = mock_model,
          U = matrix(rnorm(2 * n, sd = 0.1), n, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(rnorm(n), n, 1),
          d = matrix(rnorm(n), n, 1)
        )
      })
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  # Variance should be a 2x2 matrix
  expect_equal(dim(result$variance), c(2, 2))
  
  # Variance matrix should be symmetric
  expect_equal(result$variance, t(result$variance), tolerance = 1e-10)
  
  # Diagonal elements should be positive (variances)
  expect_true(all(diag(result$variance) > 0))
})

test_that("pool_rw computes standard errors correctly", {
  set.seed(123)
  
  mock_model <- list(coefficients = c(`(Intercept)` = 1.0, x1 = 0.5))
  class(mock_model) <- "lm"
  
  n <- 100
  with_rw_obj <- structure(
    list(
      m = 3,
      n = n,
      results = lapply(1:3, function(i) {
        list(
          model = mock_model,
          U = matrix(rnorm(2 * n, sd = 0.1), n, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(rnorm(n), n, 1),
          d = matrix(rnorm(n), n, 1)
        )
      })
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  # Standard errors should be sqrt of diagonal of variance
  expected_se <- sqrt(diag(result$variance))
  expect_equal(result$pooled$std.error, expected_se)
  
  # Standard errors should be positive
  expect_true(all(result$pooled$std.error > 0))
})

test_that("pool_rw computes test statistics correctly", {
  set.seed(123)
  
  mock_model <- list(coefficients = c(`(Intercept)` = 2.0, x1 = 1.5))
  class(mock_model) <- "lm"
  
  n <- 100
  with_rw_obj <- structure(
    list(
      m = 3,
      n = n,
      results = lapply(1:3, function(i) {
        list(
          model = mock_model,
          U = matrix(rnorm(2 * n, sd = 0.1), n, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(rnorm(n), n, 1),
          d = matrix(rnorm(n), n, 1)
        )
      })
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  # Test statistic should be estimate / std.error
  expected_stat <- result$pooled$estimate / result$pooled$std.error
  expect_equal(result$pooled$statistic, expected_stat)
})

test_that("pool_rw computes p-values correctly", {
  set.seed(123)
  
  mock_model <- list(coefficients = c(`(Intercept)` = 2.0, x1 = 1.5))
  class(mock_model) <- "lm"
  
  n <- 100
  with_rw_obj <- structure(
    list(
      m = 3,
      n = n,
      results = lapply(1:3, function(i) {
        list(
          model = mock_model,
          U = matrix(rnorm(2 * n, sd = 0.1), n, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(rnorm(n), n, 1),
          d = matrix(rnorm(n), n, 1)
        )
      })
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  # P-values should be between 0 and 1
  expect_true(all(result$pooled$p.value >= 0))
  expect_true(all(result$pooled$p.value <= 1))
  
  # P-value should be 2 * pnorm(-abs(statistic))
  expected_pval <- 2 * pnorm(-abs(result$pooled$statistic))
  expect_equal(result$pooled$p.value, expected_pval)
})

test_that("pool_rw computes confidence intervals correctly", {
  set.seed(123)
  
  mock_model <- list(coefficients = c(`(Intercept)` = 2.0, x1 = 1.5))
  class(mock_model) <- "lm"
  
  n <- 100
  with_rw_obj <- structure(
    list(
      m = 3,
      n = n,
      results = lapply(1:3, function(i) {
        list(
          model = mock_model,
          U = matrix(rnorm(2 * n, sd = 0.1), n, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(rnorm(n), n, 1),
          d = matrix(rnorm(n), n, 1)
        )
      })
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  # 95% CI should be estimate +/- 1.96 * SE
  expected_lower <- result$pooled$estimate - 1.96 * result$pooled$std.error
  expected_upper <- result$pooled$estimate + 1.96 * result$pooled$std.error
  
  expect_equal(result$pooled$conf.low, expected_lower)
  expect_equal(result$pooled$conf.high, expected_upper)
  
  # Upper bound should be greater than lower bound
  expect_true(all(result$pooled$conf.high > result$pooled$conf.low))
})

test_that("pool_rw output has correct column names", {
  mock_model <- list(coefficients = c(`(Intercept)` = 1.0, x1 = 0.5, x2 = -0.3))
  class(mock_model) <- "lm"
  
  with_rw_obj <- structure(
    list(
      m = 2,
      n = 100,
      results = lapply(1:2, function(i) {
        list(
          model = mock_model,
          U = matrix(0, 100, 3),
          tau = diag(-0.01, 3),
          S_mis_imp = matrix(0, 100, 1),
          d = matrix(0, 100, 1)
        )
      })
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  expect_named(
    result$pooled,
    c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  )
  
  expect_equal(result$pooled$term, c("(Intercept)", "x1", "x2"))
})




test_that("coef.pool_rw extracts coefficients correctly", {
  mock_model <- list(coefficients = c(`(Intercept)` = 1.5, x1 = 0.8, x2 = -0.4))
  class(mock_model) <- "lm"
  
  with_rw_obj <- structure(
    list(
      m = 2,
      n = 100,
      results = lapply(1:2, function(i) {
        list(
          model = mock_model,
          U = matrix(0, 100, 3),
          tau = diag(-0.01, 3),
          S_mis_imp = matrix(0, 100, 1),
          d = matrix(0, 100, 1)
        )
      })
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  coefs <- coef(result)
  
  expect_named(coefs, c("(Intercept)", "x1", "x2"))
  expect_equal(unname(coefs), result$pooled$estimate)
})

test_that("vcov.pool_rw returns variance matrix", {
  mock_model <- list(coefficients = c(`(Intercept)` = 1.0, x1 = 0.5))
  class(mock_model) <- "lm"
  
  n <- 100
  set.seed(123)
  with_rw_obj <- structure(
    list(
      m = 3,
      n = n,
      results = lapply(1:3, function(i) {
        list(
          model = mock_model,
          U = matrix(rnorm(2 * n, sd = 0.1), n, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(rnorm(n), n, 1),
          d = matrix(rnorm(n), n, 1)
        )
      })
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  vcov_mat <- vcov(result)
  
  expect_equal(vcov_mat, result$variance)
  expect_true(is.matrix(vcov_mat))
  expect_equal(dim(vcov_mat), c(2, 2))
})

test_that("pool_rw handles single imputation", {
  mock_model <- list(coefficients = c(`(Intercept)` = 1.0, x1 = 0.5))
  class(mock_model) <- "lm"
  
  with_rw_obj <- structure(
    list(
      m = 1,
      n = 100,
      results = list(
        list(
          model = mock_model,
          U = matrix(rnorm(200), 100, 2),
          tau = matrix(c(-0.01, 0, 0, -0.01), 2, 2),
          S_mis_imp = matrix(rnorm(100), 100, 1),
          d = matrix(rnorm(100), 100, 1)
        )
      )
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  expect_s3_class(result, "pool_rw")
  expect_equal(result$m, 1)
  expect_equal(nrow(result$pooled), 2)
})

test_that("pool_rw handles models with many predictors", {
  coef_vec <- setNames(rnorm(10), c("(Intercept)", paste0("x", 1:9)))
  mock_model <- list(coefficients = coef_vec)
  class(mock_model) <- "lm"
  
  n <- 100
  with_rw_obj <- structure(
    list(
      m = 3,
      n = n,
      results = lapply(1:3, function(i) {
        list(
          model = mock_model,
          U = matrix(rnorm(10 * n), n, 10),
          tau = diag(-0.01, 10),
          S_mis_imp = matrix(rnorm(n), n, 1),
          d = matrix(rnorm(n), n, 1)
        )
      })
    ),
    class = "with_rw"
  )
  
  result <- pool_rw(with_rw_obj)
  
  expect_equal(nrow(result$pooled), 10)
  expect_equal(dim(result$variance), c(10, 10))
})

