# rw 0.0.0.9002

- Add `pmmrw` imputation with donor recording and donor-source variance correction.
- Add PMM cross terms for a linear-model response and the documented binomial
  analysis. Other PMM analyses can supply an expected-score function.
- Preserve factor coding when computing imputation and analysis scores.
- Use separate files for imputation, analysis, donor recording, and RW pooling.
- Add GS and RW simulation scripts and task tables.

## Interface changes

- `with_rw()` returns `rw_fit`; `pool_rw()` returns `rw_pool`.
- `summary()` returns a data frame. Use `coef()` and `vcov()` for estimates
  and their covariance matrix instead of the former `$pooled` field.
- Refit saved objects from earlier versions before pooling. The unused `...`
  arguments have been removed from `with_rw()` and `pool_rw()`; the latter now
  accepts `pmm_kappa` for a supplied PMM cross term.
