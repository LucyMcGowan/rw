# Component checks adapted from the upstream helper tests.
test_that("analysis scores and Jacobians follow the fitted model and subset", {
  set.seed(71)
  data <- data.frame(x = rnorm(100))
  data$y <- 1 + data$x + rnorm(100)
  data$d <- rbinom(100, 1, plogis(data$x))
  models <- list(lm(y ~ x, data, subset = x > -1),
                 glm(d ~ x, data, subset = x > -1, family = binomial()))
  for (model in models) {
    result <- analysis_component(model, nrow(data))
    x <- model.matrix(model)
    attr(x, "assign") <- NULL
    rows <- as.integer(rownames(model.frame(model)))
    if (inherits(model, "glm")) {
      U <- x * (model.response(model.frame(model)) - fitted(model))
      tau <- -crossprod(x, x * (fitted(model) * (1 - fitted(model))))
    } else {
      U <- x * residuals(model) / summary(model)$sigma^2
      tau <- -crossprod(x) / summary(model)$sigma^2
    }
    expect_equal(unname(result$U[rows, ]), unname(U))
    expect_true(all(result$U[-rows, ] == 0))
    expect_equal(result$tau, tau)
    expect_equal(result$n_analysis, length(rows))
    expect_identical(colnames(result$U), names(coef(model)))
  }
})

test_that("parametric components use the observed-data score Jacobian", {
  data <- data.frame(x = seq(-2, 2, length.out = 8),
                     y = c(0.2, -0.8, 1.3, 0.7, -0.4, 1.8, 0.3, 1.1))
  missing <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
  for (method in c("norm", "logreg")) {
    completed <- data
    if (method == "logreg") completed$y <- c(0, 1, 0, 1, 1, 0, 1, 1)
    for (intercept in c(TRUE, FALSE)) {
      formula <- if (intercept) y ~ x else y ~ x - 1
      x <- model.matrix(formula, completed)
      attr(x, "assign") <- NULL
      p <- ncol(x)
      model <- list(setup = list(method = method), formula = formula,
                    xnames = colnames(x), beta.dot = rep(0.2, p), sigma.dot = 1.2)
      theta <- model$beta.dot
      if (method == "norm") theta <- c(theta, model$sigma.dot^2)
      score <- function(theta) {
        eta <- drop(x %*% theta[seq_len(p)])
        if (method == "logreg") return(x * (completed$y - plogis(eta)))
        variance <- theta[p + 1L]
        residual <- completed$y - eta
        cbind(x * (residual / variance),
              sigma2 = (residual^2 - variance) / (2 * variance^2))
      }
      jacobian <- vapply(seq_along(theta), function(j) {
        step <- rep(0, length(theta))
        step[j] <- 1e-6
        colSums((score(theta + step) - score(theta - step))[!missing, , drop = FALSE]) /
          (2e-6 * nrow(completed))
      }, numeric(length(theta)))
      result <- parametric_component(completed, model, "y", missing)
      expect_equal(result$S_mis_imp, score(theta) * missing)
      expect_equal(unname(result$d %*% t(jacobian)),
                   unname(-score(theta) * !missing), tolerance = 1e-7)
      expect_true(all(result$S_mis_imp[!missing, ] == 0))
      expect_true(all(result$d[missing, ] == 0))
    }
  }
})

test_that("RW assembly agrees with the component formula", {
  set.seed(72)
  n <- 8L
  for (m in c(1L, 3L)) {
    for (p in c(1L, 2L)) {
      results <- lapply(seq_len(m), function(i) {
        list(U = matrix(rnorm(n * p), n, p), tau = -diag(p) * (n + i),
             S_mis_imp = matrix(rnorm(n * 3), n, 3), d = matrix(rnorm(n * 3), n, 3))
      })
      U <- Reduce(`+`, lapply(results, `[[`, "U")) / m
      d <- Reduce(`+`, lapply(results, `[[`, "d")) / m
      kappa <- Reduce(`+`, lapply(results, function(x) crossprod(x$U, x$S_mis_imp))) / (m * n)
      alpha <- Reduce(`+`, lapply(results, function(x) crossprod(x$d))) / (m * n)
      tau <- Reduce(`+`, lapply(results, `[[`, "tau")) / (m * n)
      cross <- kappa %*% crossprod(d, U) / n
      middle <- crossprod(U) / n + kappa %*% alpha %*% t(kappa) + cross + t(cross)
      expected <- solve(tau) %*% middle %*% t(solve(tau)) / n
      variance <- rw_variance(list(results = results, m = m, n = n))
      expect_equal(variance, expected)
      expect_equal(variance, t(variance))
      expect_equal(dim(variance), c(p, p))
    }
  }
})

test_that("pooling and inference preserve coefficients and use model-specific quantiles", {
  set.seed(73)
  data <- data.frame(x = rnorm(80), y = rnorm(80))
  data$d <- rbinom(80, 1, plogis(data$x))
  for (models in list(list(lm(y ~ 1, data)),
                      list(lm(y ~ x, data[1:40, ]), lm(y ~ x, data[41:80, ])),
                      list(glm(d ~ x, data, family = binomial())))) {
    results <- lapply(models, function(model) {
      c(list(model = model), analysis_component(model, nrow(data)),
        list(S_mis_imp = matrix(0, nrow(data), 1), d = matrix(0, nrow(data), 1)))
    })
    fit <- structure(list(results = results, m = length(models), n = nrow(data),
                          mids = list(method = character(), models = list())), class = "rw_fit")
    pooled <- pool_rw(fit)
    expected <- Reduce(`+`, lapply(models, coef)) / length(models)
    expect_equal(coef(pooled), expected)
    table <- summary(pooled)
    se <- sqrt(diag(vcov(pooled)))
    statistic <- expected / se
    if (inherits(models[[1L]], "glm")) {
      critical <- qnorm(0.975)
      p_value <- 2 * pnorm(-abs(statistic))
    } else {
      critical <- qt(0.975, df.residual(models[[1L]]))
      p_value <- 2 * pt(-abs(statistic), df.residual(models[[1L]]))
    }
    expect_named(table, c("term", "estimate", "std.error", "statistic", "p.value",
                          "conf.low", "conf.high"))
    expect_equal(table$term, names(expected))
    expect_equal(table$statistic, unname(statistic))
    expect_equal(table$p.value, unname(p_value))
    expect_equal(table$conf.low, unname(expected - critical * se))
    expect_equal(table$conf.high, unname(expected + critical * se))
  }
})
