# Fit analysis model with Robins-Wang variance calculation

Fit analysis model with Robins-Wang variance calculation

## Usage

``` r
with_rw(data, expr, ...)
```

## Arguments

- data:

  A mids object from mice with parametric imputation methods

- expr:

  Expression to evaluate on each imputed dataset (typically a model
  formula)

- ...:

  Additional arguments passed to the analysis function

## Value

A with_rw object containing fitted models and RW components

## Details

IMPORTANT: This function requires parametric imputation methods that
produce score functions and information matrices. Supported methods
include:

- norm (for continuous variables)

- logreg (for binary variables)

## Examples

``` r
if (FALSE) { # \dontrun{
library(mice)
imp <- mice(nhanes, method = "norm", m = 5, tasks = "train", print = FALSE)
fit <- with_rw(imp, lm(bmi ~ age + hyp))
pool_rw(fit)
} # }
```
