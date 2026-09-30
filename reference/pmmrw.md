# PMM imputation and donor recording

Use the PMM donor draw implemented in `mice` and record donor IDs and
fitted matching models.

## Usage

``` r
mice.impute.pmmrw(y, ry, x, wy = !ry, task = "impute", model = NULL,
  exclude = NULL, ridge = 1e-05, matchtype = 1L, donors = 5L,
  use.matcher = FALSE, mlocal = 1L, ...)
extract_donor_id(data, variable)
```

## Arguments

- y, ry, x, wy:

  Data and row indicators supplied by `mice`.

- task, model:

  Recording arguments supplied by `mice` with `tasks = "train"`.

- exclude, matchtype, use.matcher, mlocal:

  Keep `exclude = NULL`, `matchtype = 1`, `use.matcher = FALSE` and
  `mlocal = 1`.

- ridge:

  Ridge parameter used by `mice`.

- donors:

  Donor pool size, default 5.

- ...:

  Additional arguments passed to `mice`'s normal draw.

- data:

  A `mids` object containing a `pmmrw` imputation.

- variable:

  Name of the PMM-imputed variable.

## Value

`mice.impute.pmmrw()` returns values copied from the selected donors.
`extract_donor_id()` returns an `n` by `m` matrix: rows follow the
original data and columns follow imputations. Imputed rows contain the
selected donor's original row number; observed rows contain `NA`.

## Details

Set `method = "pmmrw"` for a numeric variable in
`mice(..., tasks = "train")`.

With the supported PMM settings described above, using the same random
seed and donor pool produces the same imputed values as
[`mice.impute.pmm()`](https://amices.org/mice/reference/mice.impute.pmm.html).

## See also

[`mice.impute.pmm`](https://amices.org/mice/reference/mice.impute.pmm.html)
for PMM arguments.
