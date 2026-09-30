# PMM cross terms

Combine expected analysis scores with derivatives of donor-selection
probabilities. Matching-model predictors must be fully observed.

## Usage

``` r
pmm_kappa(object, variable, expected_score, ...)
pmm_kappa_binomial(object, variable, threshold = -Inf, quadrature_order = 24L)
```

## Arguments

- object:

  An `rw_fit` object from
  [`with_rw`](https://lucymcgowan.github.io/rw/reference/rw.md).

- variable:

  Name of the numeric PMM-imputed variable.

- expected_score:

  Function returning expected analysis-score matrices; see below.

- ...:

  Additional arguments passed to `expected_score`.

- threshold:

  Lower analysis threshold. Match the subset `A > threshold`; the
  default `-Inf` uses all rows.

- quadrature_order:

  Number of quadrature points, default 24.

## Value

The PMM columns of the RW cross-term matrix. Pass this matrix to
[`pool_rw`](https://lucymcgowan.github.io/rw/reference/rw.md) as
`pmm_kappa`.

## Expected scores

`pmm_kappa()` calls
`expected_score(object, variable, p, missing, donor_mean, ...)`. Here
`p` is the imputation index, `missing` is a length-`n` logical indicator
and `donor_mean` contains fitted donor predictions.

Return one matrix per analysis coefficient, in coefficient order. Rows
are observed donors and columns are missing recipients, both in original
data order. Each entry is the expected score for that donor-recipient
pair, accounting for the analysis subset and downstream imputation.

The cross term averages scores times probability derivatives over
imputations and divides by the full sample size. Only the selection
probabilities are differentiated, not the expected scores.

## Binomial analysis

`pmm_kappa_binomial()` supports an unweighted
`glm(D ~ A, family = binomial())` with the default logit link,
optionally using `A > threshold`. `A` uses `pmmrw` and `D` uses
`logreg`. The logistic imputation model must include `A` as an
untransformed additive predictor; its other predictors must be fully
observed.

The score uses observed responses where available and the recorded
logistic model otherwise. Quadrature integrates over the fitted normal
donor distribution; it does not change the donor pool or imputed values.
