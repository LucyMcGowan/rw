# Setup: Create test data and mice objects for testing
library(mice)

# Helper to create a dataset with relationships between variables
create_test_data <- function(n = 100, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 0.5 + 1.2 * x1 - 0.8 * x2 + rnorm(n, sd = 0.5)
  
  data.frame(
    x1 = x1,
    x2 = x2,
    y = y
  )
}

# Helper to introduce missing data
add_missing_data <- function(data, vars = "y", prop = 0.2) {
  for (var in vars) {
    n_miss <- round(nrow(data) * prop)
    miss_idx <- sample(nrow(data), n_miss)
    data[[var]][miss_idx] <- NA
  }
  data
}

# Create mice object with parametric methods
create_mice_object <- function(data, method = "norm", m = 3, seed = 1) {
  set.seed(seed)
  suppressWarnings(
    mice(data, method = method, m = m, maxit = 5, 
         printFlag = FALSE, seed = seed, tasks = "train")
  )
}

test_that("with_rw requires mids object", {
  expect_error(
    with_rw(data.frame(x = 1:5), lm(y ~ x))
  )
  
  expect_error(
    with_rw(list(m = 5), lm(y ~ x))
  )
})

test_that("with_rw requires imputation models", {
  data <- create_test_data()
  data <- add_missing_data(data, "y")
  
  # Create mice object without saving models
  imp <- suppressWarnings(
    mice(data, method = "norm", m = 3, maxit = 5, printFlag = FALSE)
  )
  # imp will not have models since tasks != "train"
  
  expect_error(
    with_rw(imp, lm(y ~ x1 + x2)),
    "must contain imputation models"
  )
})

test_that("with_rw checks for supported imputation methods - pmm", {
  data <- create_test_data()
  data <- add_missing_data(data, "y")
  
  imp <- suppressWarnings(
    mice(data, method = "pmm", m = 3, maxit = 5, 
         printFlag = FALSE, tasks = "train")
  )
  
  expect_error(
    with_rw(imp, lm(y ~ x1 + x2)),
    "Unsupported imputation method"
  )
  
  expect_error(
    with_rw(imp, lm(y ~ x1 + x2)),
    "pmm"
  )
})

test_that("with_rw checks for supported imputation methods - cart", {
  data <- create_test_data()
  data <- add_missing_data(data, "y")
  
  imp <- suppressWarnings(
    mice(data, method = "cart", m = 3, maxit = 5, 
         printFlag = FALSE, tasks = "train")
  )
  
  expect_error(
    with_rw(imp, lm(y ~ x1 + x2)),
    "cart"
  )
})

test_that("with_rw accepts norm method", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y", prop = 0.2)
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  
  expect_no_error({
    result <- with_rw(imp, lm(y ~ x1 + x2))
  })
})

test_that("with_rw returns with_rw object", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 3)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_s3_class(result, "with_rw")
})

test_that("with_rw object has correct structure", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 3)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_named(result, c("call", "results", "m", "n", "expr", "mids"))
  expect_equal(result$m, 3)
  expect_equal(result$n, 50)
  expect_length(result$results, 3)
})

test_that("with_rw processes each imputation", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 5)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_length(result$results, 5)
  
  # Each result should have the required components
  for (i in 1:5) {
    expect_named(result$results[[i]], 
                 c("model", "U", "tau", "S_mis_imp", "d", "n_analysis"))
  }
})

test_that("with_rw stores correct number of observations", {
  data <- create_test_data(n = 80)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 3)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_equal(result$n, 80)
})

test_that("with_rw fits analysis models correctly", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 3)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  # Check that models were fit
  for (i in 1:3) {
    expect_s3_class(result$results[[i]]$model, "lm")
    expect_named(coef(result$results[[i]]$model), 
                 c("(Intercept)", "x1", "x2"))
  }
})

test_that("with_rw handles lm analysis models", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 3)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  # Check U matrix dimensions
  expect_equal(nrow(result$results[[1]]$U), 50)
  expect_equal(ncol(result$results[[1]]$U), 3)  # intercept + x1 + x2
  
  # Check tau matrix dimensions
  expect_equal(dim(result$results[[1]]$tau), c(3, 3))
})

test_that("with_rw handles glm analysis models", {
  set.seed(1)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  prob <- plogis(0.5 + 0.8 * x1 - 0.6 * x2)
  y <- rbinom(n, 1, prob)
  data <- data.frame(x1 = x1, x2 = x2,  y = factor(y, levels = c(0, 1)))
  data <- add_missing_data(data, "y")
  
  methods <- c(x1 = "", x2 = "", y = "logreg") 

  imp <- create_mice_object(data, method = methods, m = 3)
  result <- with_rw(imp, glm(y ~ x1 + x2, family = binomial))
  
  # Check that GLMs were fit
  for (i in 1:3) {
    expect_s3_class(result$results[[i]]$model, "glm")
  }
  
  # Check U matrix dimensions
  expect_equal(nrow(result$results[[1]]$U), 50)
  expect_equal(ncol(result$results[[1]]$U), 3)
})

test_that("with_rw requires lm or glm analysis model", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  
  expect_error(
    with_rw(imp, mean(y)),
    "Analysis model must be.*lm.*or.*glm"
  )
})

test_that("with_rw extends U to full dataset", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  # U should have n rows (full dataset)
  expect_equal(nrow(result$results[[1]]$U), 50)
})


test_that("with_rw computes S_mis_imp correctly", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  # S_mis_imp should have n rows
  expect_equal(nrow(result$results[[1]]$S_mis_imp), 50)
  
  # Should be numeric matrix
  expect_true(is.matrix(result$results[[1]]$S_mis_imp))
  expect_true(is.numeric(result$results[[1]]$S_mis_imp))
})

test_that("with_rw computes d matrix correctly", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  # d should have n rows
  expect_equal(nrow(result$results[[1]]$d), 50)
  
  # Should be numeric matrix
  expect_true(is.matrix(result$results[[1]]$d))
  expect_true(is.numeric(result$results[[1]]$d))
})

test_that("with_rw stores n_analysis", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y", prop = 0.2)
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  # n_analysis should be less than or equal to n
  expect_true(result$results[[1]]$n_analysis <= 50)
  expect_true(result$results[[1]]$n_analysis > 0)
})

test_that("with_rw works with single imputed variable", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 3)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_s3_class(result, "with_rw")
  expect_length(result$results, 3)
})

test_that("with_rw works with multiple imputed variables", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, c("x1", "y"), prop = 0.2)
  
  imp <- create_mice_object(data, method = "norm", m = 3)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_s3_class(result, "with_rw")
  expect_length(result$results, 3)
})

test_that("with_rw preserves call", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_true(is.call(result$call))
})

test_that("with_rw preserves expression", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_true(!is.null(result$expr))
})

test_that("with_rw stores mids object reference", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_s3_class(result$mids, "mids")
  expect_identical(result$mids, imp)
})

test_that("with_rw handles different sample sizes", {
  # Small sample
  data_small <- create_test_data(n = 30)
  data_small <- add_missing_data(data_small, "y")
  imp_small <- create_mice_object(data_small, method = "norm", m = 2)
  result_small <- with_rw(imp_small, lm(y ~ x1 + x2))
  
  expect_equal(result_small$n, 30)
  
  # Large sample
  data_large <- create_test_data(n = 200)
  data_large <- add_missing_data(data_large, "y")
  imp_large <- create_mice_object(data_large, method = "norm", m = 2)
  result_large <- with_rw(imp_large, lm(y ~ x1 + x2))
  
  expect_equal(result_large$n, 200)
})

test_that("with_rw handles different numbers of imputations", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  # m = 2
  imp2 <- create_mice_object(data, method = "norm", m = 2)
  result2 <- with_rw(imp2, lm(y ~ x1 + x2))
  expect_equal(result2$m, 2)
  expect_length(result2$results, 2)
  
  # m = 10
  imp10 <- create_mice_object(data, method = "norm", m = 10)
  result10 <- with_rw(imp10, lm(y ~ x1 + x2))
  expect_equal(result10$m, 10)
  expect_length(result10$results, 10)
})

test_that("with_rw handles models with intercept", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2))
  
  expect_true("(Intercept)" %in% names(coef(result$results[[1]]$model)))
})

test_that("with_rw handles models without intercept", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2 - 1))
  
  expect_false("(Intercept)" %in% names(coef(result$results[[1]]$model)))
})

test_that("with_rw handles single predictor", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1))
  
  expect_equal(ncol(result$results[[1]]$U), 2)  # intercept + x1
})

test_that("with_rw handles many predictors", {
  set.seed(1)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  x5 <- rnorm(n)
  y <- 0.5 + 0.3*x1 + 0.2*x2 - 0.4*x3 + 0.1*x4 - 0.3*x5 + rnorm(n, sd = 0.5)
  
  data <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, y = y)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 2)
  result <- with_rw(imp, lm(y ~ x1 + x2 + x3 + x4 + x5))
  
  expect_equal(ncol(result$results[[1]]$U), 6)  # intercept + 5 predictors
})

test_that("with_rw validates all imputation methods", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, c("x1", "y"))
  
  # Create imputation with mixed methods - specify pmm for both
  meth <- c("pmm", "pmm", "")
  names(meth) <- c("x1", "y", "x2")
  
  imp <- suppressWarnings(
    mice(data, method = meth, m = 2, maxit = 5, 
         printFlag = FALSE, seed = 123, tasks = "train")
  )
  
  expect_error(
    with_rw(imp, lm(y ~ x1 + x2)),
    "pmm"
  )
})

test_that("with_rw error message shows unsupported methods", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- suppressWarnings(
    mice(data, method = c("", "", "rf"), m = 2, 
         printFlag = FALSE, tasks = "train")
  )
  
  expect_error(
    with_rw(imp, lm(y ~ x1 + x2)),
    "rf"
  )
})

test_that("with_rw error message shows supported methods", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- suppressWarnings(
    mice(data, method = "cart", m = 2, maxit = 5, 
         printFlag = FALSE, tasks = "train")
  )
  
  expect_error(
    with_rw(imp, lm(y ~ x1 + x2)),
    "norm"
  )
  
  expect_error(
    with_rw(imp, lm(y ~ x1 + x2)),
    "logreg"
  )
})

test_that("with_rw errors for norm.nob method", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm.nob", m = 2)
  
  expect_error({
    result <- with_rw(imp, lm(y ~ x1 + x2))
  })
})

test_that("with_rw errors norm.boot method", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm.boot", m = 2)
  
  expect_error({
    result <- with_rw(imp, lm(y ~ x1 + x2))
  })
})

test_that("with_rw errors norm.predict method", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm.predict", m = 2)
  
  expect_error({
    result <- with_rw(imp, lm(y ~ x1 + x2))
  })
})



test_that("with_rw works end-to-end with pool_rw", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 3)
  fit <- with_rw(imp, lm(y ~ x1 + x2))
  pooled <- pool_rw(fit)
  
  expect_s3_class(pooled, "pool_rw")
  expect_equal(nrow(pooled$pooled), 3)  # intercept + x1 + x2
})

test_that("with_rw produces consistent results across runs with same seed", {
  data <- create_test_data(n = 50)
  data <- add_missing_data(data, "y")
  
  imp <- create_mice_object(data, method = "norm", m = 3, seed = 456)
  result1 <- with_rw(imp, lm(y ~ x1 + x2))
  
  imp2 <- create_mice_object(data, method = "norm", m = 3, seed = 456)
  result2 <- with_rw(imp2, lm(y ~ x1 + x2))
  
  # Coefficients should be identical
  coef1 <- sapply(result1$results, function(r) coef(r$model))
  coef2 <- sapply(result2$results, function(r) coef(r$model))
  
  expect_equal(coef1, coef2)
})

