# Pool results with Robins-Wang variance

Pool results with Robins-Wang variance

## Usage

``` r
pool_rw(object, ...)
```

## Arguments

- object:

  A with_rw object from with_rw()

- ...:

  Additional arguments (currently unused)

## Value

A pool_rw object with pooled estimates and Robins-Wang SEs

## Examples

``` r
if (FALSE) { # \dontrun{
library(mice)
imp <- mice(nhanes, method = "norm", m = 5, tasks = "train", print = FALSE)
fit <- with_rw(imp, lm(bmi ~ age + hyp))
pooled <- pool_rw(fit)
summary(pooled)
} # }
```
