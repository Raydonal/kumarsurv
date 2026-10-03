# Simulate from the Gaussian-copula marginal regression

Used for Monte Carlo validation and for the parametric bootstrap behind
the paper's goodness-of-fit tests and prediction bands.

## Usage

``` r
sim_copula(
  X,
  Z,
  beta,
  gamma,
  ar = numeric(0),
  ma = numeric(0),
  family = "kuma",
  tau = 0.5
)
```

## Arguments

- X:

  Location-submodel design matrix.

- Z:

  Shape-submodel design matrix.

- beta, gamma:

  Numeric vectors of location/shape-submodel coefficients.

- ar, ma:

  Numeric vectors of ARMA coefficients (empty for independence).

- family:

  `"kuma"` (default) or `"beta"`.

- tau:

  Quantile level for the Kumaraswamy margin (ignored if
  `family = "beta"`).

## Value

Numeric vector of length `nrow(X)`, a simulated series in `(0, 1)`.
