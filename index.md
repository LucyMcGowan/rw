# rw

The purpose of the rw package is to compute Robins-Wang variance
estimates for multiply imputed data analyses.

## Installation

This package requires the development version of mice with the `tasks`
argument:

``` r
remotes::install_github("amices/mice")
```

Then you can install the development version of rw like so:

``` r
remotes::install_github("LucyMcGowan/rw")
```

## Example

The first step is to impute your data using mice. Note that currently
only `method = "norm"` and `method = "logreg"` are supported. When using
the `mice` function, you must set the parameter `task = "train"` in
order to pass the necessary parts to our subsequent functions. For
example, see below. We will use the `nhanes` data from the mice package
and perform multiple imputation using the `norm` method.

``` r
library(rw)
library(mice)
imp <- mice(nhanes, method = "norm", m = 5, tasks = "train", print = FALSE)
```

Now, suppose we want to fit a model predicting `bmi` from `age` and
`hyp`. The `with_rw` function allows this, where the first argument is
the imputation object from the `mice` function above and the second is
the model expression. Note that currently outcome models must be
Gaussian or Binomial.

``` r
fit <- with_rw(imp, lm(bmi ~ age + hyp))
```

Finally, we can pool the results from the fit object and calculate the
Robins-Wang variance using the `pool_rw` function

``` r
pool_rw(fit)
#> 
#> ── Robins-Wang Pooled Results ──────────────────────────────────────────────────
#> Number of imputations: 5
#> Sample size: 25
#> 
#>          term estimate std.error statistic p.value conf.low conf.high
#> 1 (Intercept)   28.617     51.88   0.55161  0.5812   -73.06    130.30
#> 2         age   -3.058     23.33  -0.13105  0.8957   -48.79     42.67
#> 3         hyp    3.204     34.36   0.09325  0.9257   -64.13     70.54
```

Let’s compare this result to using Rubin’s rules.

``` r
fit_rr <- with(imp, lm(bmi ~ age + hyp))
pool(fit_rr) |>
  summary()
#>          term  estimate std.error statistic        df      p.value
#> 1 (Intercept) 28.616514  3.403986  8.406766 11.558784 2.911276e-06
#> 2         age -3.057764  1.436188 -2.129083  8.956052 6.226165e-02
#> 3         hyp  3.203510  2.961103  1.081864  7.473962 3.129791e-01
```

## Giganti & Shepherd Example

Below is an example to replicate the Giganti & Shepherd (2020) results.

``` r
set.seed(1)
meth <- make.method(giganti_data)
meth[] <- ""
meth[c("A", "D")] <- c("norm", "logreg")

pred <- make.predictorMatrix(giganti_data)
pred[,] <- 0
pred["A", c("X1", "X2", "A.star", "D.star")] <- 1
pred["D", c("X1", "X2", "A.star", "D.star", "A")] <- 1

imp <- mice(giganti_data, m = 10, method = meth, predictorMatrix = pred,
            print = FALSE, tasks = "train")

fit_rw <- with_rw(imp, glm(D ~ A, subset = A > 2, family = binomial()))

pooled <- pool_rw(fit_rw)
pooled
#> 
#> ── Robins-Wang Pooled Results ──────────────────────────────────────────────────
#> Number of imputations: 10
#> Sample size: 4000
#> 
#>          term estimate std.error statistic   p.value conf.low conf.high
#> 1 (Intercept)  -3.4864   0.23876    -14.60 2.721e-48  -3.9543    -3.018
#> 2           A   0.8764   0.08217     10.67 1.469e-26   0.7154     1.038
```
