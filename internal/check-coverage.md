# Check Coverage


``` r
library(mice)
library(rw)
library(dplyr)
library(purrr)


simulate_once <- function(i, n, m, true_beta, sigma, miss_prop) {
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- true_beta[1] + true_beta[2] * X1 + true_beta[3] * X2 + rnorm(n, sd = sigma)
  
  dat <- data.frame(Y = Y, X1 = X1, X2 = X2)
  miss_idx <- sample(1:n, size = floor(miss_prop * n))
  miss_idx_1 <- sample(1:n, size = floor(miss_prop * n))
  dat$Y[miss_idx] <- NA
  dat$X1[miss_idx_1] <- NA
  
  imp <- mice(dat, method = "norm", m = m, tasks = "train", printFlag = FALSE)
  
  rw_result <- tryCatch({
    fit_rw <- with_rw(imp, lm(Y ~ X1 + X2))
    pooled_rw <- pool_rw(fit_rw)
    
    data.frame(
      sim = i,
      method = "RW",
      term = pooled_rw$pooled$term,
      estimate = pooled_rw$pooled$estimate,
      se = pooled_rw$pooled$std.error
    )
  }, error = function(e) NULL)
  
  rubin_result <- tryCatch({
    fit_mice <- with(imp, lm(Y ~ X1 + X2))
    pooled_mice <- pool(fit_mice)
    summ <- summary(pooled_mice)
    
    data.frame(
      sim = i,
      method = "Rubin",
      term = summ$term,
      estimate = summ$estimate,
      se = summ$std.error
    )
  }, error = function(e) NULL)
  
  bind_rows(rw_result, rubin_result)
}
```

``` r
n_sims <- 500

set.seed(1)

true_params <- data.frame(
  term = c("(Intercept)", "X1", "X2"),
  true_value = c(2, 0, -1)
)
n <- 1000
crit <- qt(0.975, n - ncol(true_params))
results <- seq_len(n_sims) |>
  map(simulate_once, n = n, m = 20, 
      true_beta = true_params$true_value, 
      sigma = 5, miss_prop = 0.3) |>
  list_rbind()


results |>
  left_join(true_params, by = "term") |>
  group_by(method, term) |>
  summarise(
    empirical_sd = sd(estimate),
    estiamted_sd = mean(se),
    coverage_95 = mean(abs(estimate - true_value) <= crit * se),
    .groups = "drop"
  )
```

    # A tibble: 6 × 5
      method term        empirical_sd estiamted_sd coverage_95
      <chr>  <chr>              <dbl>        <dbl>       <dbl>
    1 RW     (Intercept)        0.197        0.189       0.942
    2 RW     X1                 0.231        0.211       0.94 
    3 RW     X2                 0.183        0.189       0.968
    4 Rubin  (Intercept)        0.197        0.190       0.938
    5 Rubin  X1                 0.231        0.228       0.95 
    6 Rubin  X2                 0.183        0.190       0.964
