# Giganti and Shepherd simulation data

A simulated two-phase dataset based on Giganti and Shepherd (2020), not
participant data from their study. Validated values of `A` and `D` are
available for 1,000 of the 4,000 observations.

## Usage

``` r
giganti_data
```

## Format

A data frame with 4,000 rows and 8 variables:

- ID:

  Observation identifier.

- X1, X2:

  Fully observed continuous covariates.

- A.star:

  Fully observed error-prone measure of the continuous exposure.

- D.star:

  Fully observed error-prone binary outcome, coded 0/1.

- A:

  Validated continuous exposure; missing outside the validation sample.

- D:

  Validated binary outcome, a factor with levels 0 and 1; missing
  outside the validation sample.

- Intercept:

  A constant column equal to 1.

## References

Giganti MJ, Shepherd BE (2020). Multiple-Imputation Variance Estimation
in Studies With Missing or Misclassified Inclusion Criteria. *American
Journal of Epidemiology*, 189(12), 1628–1632.
[doi:10.1093/aje/kwaa153](https://doi.org/10.1093/aje/kwaa153) .
