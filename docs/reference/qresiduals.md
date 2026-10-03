# Predictive quantile residuals

One-step-ahead predictive quantile residuals (Masarotto and Varin, 2012,
2017) from a model fitted by
[`fit_copula()`](https://raydonal.github.io/kumarsurv/reference/fit_copula.md):
the standardized innovations evaluated at the maximum-likelihood
estimates. Under correct specification and known parameters, these are
i.i.d. standard normal; see the paper's Proposition 1 for the
finite-sample statement.

## Usage

``` r
qresiduals(fit)
```

## Arguments

- fit:

  A list as returned by
  [`fit_copula()`](https://raydonal.github.io/kumarsurv/reference/fit_copula.md).

## Value

A list with components `r` (the residuals), `mu`, `shape` (the fitted
marginal parameters) and `eps` (the untransformed latent scale).

## See also

[`fit_copula()`](https://raydonal.github.io/kumarsurv/reference/fit_copula.md)
