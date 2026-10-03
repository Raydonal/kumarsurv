# Fit the Gaussian-copula marginal regression by maximum likelihood

Fits the model of
[`nll_copula()`](https://raydonal.github.io/kumarsurv/reference/nll_copula.md)
by BFGS (Nelder-Mead is unstable on this likelihood surface; see the
paper's pitfall 6), starting from the independence-model fit, then
computes standard errors from the numerical Hessian with a relative step
of `1e-3` (the `numDeriv` default, `d = 0.1`, can push an autoregressive
coefficient outside the stationary region near the boundary; see the
paper's pitfall 16).

## Usage

``` r
fit_copula(y, X, Z, p = 0, q = 0, family = "kuma", tau = 0.5)
```

## Arguments

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

A list with components `par` (named numeric vector), `se` (standard
errors, `NA` if the Hessian is not invertible), `loglik`, `npar`, `aic`,
`bic`, `family`, `p`, `q`, `tau`, the inputs `X`, `Z`, `y`, and
`convergence` (as returned by
[`stats::optim()`](https://rdrr.io/r/stats/optim.html)).

## Details

**At `tau` other than the median, prefer
[`fit_copula_multi_tau()`](https://raydonal.github.io/kumarsurv/reference/fit_copula_multi_tau.md)**:
the likelihood surface can have more than one interior stationary point
at a given `tau`, with different values, and a bare single-start
`fit_copula()` call is not guaranteed to land on the one with the higher
likelihood (see the paper's pitfall 1).

## See also

[`fit_copula_multi_tau()`](https://raydonal.github.io/kumarsurv/reference/fit_copula_multi_tau.md),
[`qresiduals()`](https://raydonal.github.io/kumarsurv/reference/qresiduals.md),
[`sim_copula()`](https://raydonal.github.io/kumarsurv/reference/sim_copula.md)
