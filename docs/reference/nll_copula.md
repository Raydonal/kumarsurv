# Negative log-likelihood of the Gaussian-copula marginal regression

The likelihood underlying
[`fit_copula()`](https://raydonal.github.io/kumarsurv/reference/fit_copula.md):
continuous Kumaraswamy or beta margins, serial dependence induced by a
Gaussian copula with ARMA(`p`,`q`) structure on the latent
(probit-transformed) scale. Evaluated in closed form via the innovations
decomposition (no explicit construction or inversion of the `n x n`
correlation matrix).

## Usage

``` r
nll_copula(par, y, X, Z, p, q, family = "kuma", tau = 0.5)
```

## Arguments

- par:

  Numeric vector of parameters, in the order: location-submodel
  coefficients (length `ncol(X)`), shape-submodel coefficients (length
  `ncol(Z)`), AR coefficients (length `p`), MA coefficients (length
  `q`).

- y:

  Numeric response vector in `(0, 1)`.

- X:

  Location-submodel design matrix.

- Z:

  Shape-submodel design matrix.

- p, q:

  Non-negative integers, the ARMA orders.

- family:

  `"kuma"` (default) or `"beta"`.

- tau:

  Quantile level for the Kumaraswamy margin (ignored if
  `family = "beta"`).

## Value

The negative log-likelihood, a single number (`1e10` as a penalty for
parameter values outside the admissible region, e.g. a non-stationary
ARMA polynomial).

## See also

[`fit_copula()`](https://raydonal.github.io/kumarsurv/reference/fit_copula.md)
