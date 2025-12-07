# Check Coverage


``` r
library(mice)
```


    Attaching package: 'mice'

    The following object is masked from 'package:stats':

        filter

    The following objects are masked from 'package:base':

        cbind, rbind

``` r
library(rw)
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
library(purrr)


simulate_once <- function(i, n, m, true_beta, sigma, miss_prop) {
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- true_beta[1] + true_beta[2] * X1 + true_beta[3] * X2 + rnorm(n, sd = sigma)
  
  dat <- data.frame(Y = Y, X1 = X1, X2 = X2)
  miss_idx <- sample(1:n, size = floor(miss_prop * n))
  dat$Y[miss_idx] <- NA
  
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
n_sims <- 1000

set.seed(1)

true_params <- data.frame(
  term = c("(Intercept)", "X1", "X2"),
  true_value = c(2, 0, -1)
)
results <- seq_len(n_sims) |>
  map(simulate_once, n = 1000, m = 10, 
      true_beta = true_params$true_value, 
      sigma = 2, miss_prop = 0.1) |>
  list_rbind()


results |>
  left_join(true_params, by = "term") |>
  group_by(method, term) |>
  summarise(
    empirical_sd = sd(estimate),
    estiamted_sd = mean(se),
    coverage_95 = mean(abs(estimate - true_value) <= 1.96 * se),
    .groups = "drop"
  )
```

    # A tibble: 6 × 5
      method term        empirical_sd estiamted_sd coverage_95
      <chr>  <chr>              <dbl>        <dbl>       <dbl>
    1 RW     (Intercept)       0.0679       0.0683       0.951
    2 RW     X1                0.0686       0.0748       0.968
    3 RW     X2                0.0701       0.0711       0.95 
    4 Rubin  (Intercept)       0.0679       0.0669       0.949
    5 Rubin  X1                0.0686       0.0670       0.942
    6 Rubin  X2                0.0701       0.0669       0.937
