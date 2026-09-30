# Fit and pool an RW analysis

Fit an analysis to each completed dataset and combine the estimates
using Robins-Wang variance estimation.

## Usage

``` r
with_rw(data, expr)
pool_rw(object, pmm_kappa = NULL)
```

## Arguments

- data:

  A `mids` object from `mice(..., tasks = "train")` using `norm`,
  `logreg`, or `pmmrw`.

- expr:

  An unweighted `lm` or binomial `glm` expression. Use a numeric
  response for `lm`; for `glm`, use the default logit link and a 0/1 or
  two-level factor response.

- object:

  An `rw_fit` object from `with_rw()`.

- pmm_kappa:

  Optional PMM cross-term matrix.

## Value

`with_rw()` returns an `rw_fit` with the fitted models, RW components
and original imputations. `pool_rw()` returns an `rw_pool` with pooled
estimates and variance. Use
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) to extract results.

## Details

`with_rw()` retains completed values and variable types. Rows excluded
by the analysis subset receive zero analysis scores.

Parametric cross terms are calculated from the stored analysis and
imputation scores. PMM pooling supports one numeric PMM variable with
fully observed matching-model predictors and accounts for repeated donor
use.

Automatic PMM pooling requires the untransformed PMM variable as the
[`lm()`](https://rdrr.io/r/stats/lm.html) response. Analysis predictors
and any subset must depend only on fully observed variables. Other
analyses need a matrix from
[`pmm_kappa`](https://lucymcgowan.github.io/rw/reference/pmm.md) or
[`pmm_kappa_binomial`](https://lucymcgowan.github.io/rw/reference/pmm.md).

[`summary()`](https://rdrr.io/r/base/summary.html) reports standard
errors, tests and 95% confidence intervals. It uses the first fitted
model's residual degrees of freedom for
[`lm()`](https://rdrr.io/r/stats/lm.html) and the normal distribution
for binomial [`glm()`](https://rdrr.io/r/stats/glm.html).
