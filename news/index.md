# Changelog

## rw 0.1.0

Initial CRAN submission.

## rw 0.0.0.9002

- Add `pmmrw` imputation with donor recording and donor-source variance
  correction.
- Add PMM cross terms for a linear-model response and the documented
  binomial analysis. Other PMM analyses can supply an expected-score
  function.
- Preserve factor coding when computing imputation and analysis scores.
- Use separate files for imputation, analysis, donor recording, and RW
  pooling.
- Add GS and RW simulation scripts and task tables.

### Interface changes

- [`with_rw()`](https://lucymcgowan.github.io/rw/reference/rw.md)
  returns `rw_fit`;
  [`pool_rw()`](https://lucymcgowan.github.io/rw/reference/rw.md)
  returns `rw_pool`.
- [`summary()`](https://rdrr.io/r/base/summary.html) returns a data
  frame. Use [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) for estimates and their
  covariance matrix instead of the former `$pooled` field.
- Refit saved objects from earlier versions before pooling. The unused
  `...` arguments have been removed from
  [`with_rw()`](https://lucymcgowan.github.io/rw/reference/rw.md) and
  [`pool_rw()`](https://lucymcgowan.github.io/rw/reference/rw.md); the
  latter now accepts `pmm_kappa` for a supplied PMM cross term.
