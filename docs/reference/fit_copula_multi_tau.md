# Fit the model at several quantile levels, with a monotonicity/invariance check

Fits
[`fit_copula()`](https://raydonal.github.io/kumarsurv/reference/fit_copula.md)
at each level in `taus` and runs two systematic diagnostics the paper
relies on to avoid a known failure mode (pitfall 1: a single-start fit
at `tau = 0.90` can land on a spurious stationary point whose fitted
quantile sits above every observation):

## Usage

``` r
fit_copula_multi_tau(
  y,
  X,
  Z,
  p = 0,
  q = 0,
  family = "kuma",
  taus = c(0.5, 0.75, 0.9, 0.95, 0.99),
  tol_ll = 1,
  verbose = TRUE
)
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

- taus:

  Numeric vector of quantile levels to fit.

- tol_ll:

  Absolute log-likelihood tolerance for the invariance check (default
  `1`; fits whose log-likelihood differs from the median by more than
  this are flagged).

- verbose:

  Logical; print the diagnostic table.

## Value

A list named by `tau`, each element with components `par`, `loglik`,
`fonte` (`"padrao"` or `"reinicializado"`: standard fit or restarted)
and `aprovado` (logical, passed both checks), with a `"diagnostico"`
attribute (a data frame) summarizing the checks for every level.

## Details

1.  **Log-likelihood invariance.** The `tau`-reparametrization describes
    the *same* distribution family for every `tau`; only the
    parametrization changes, not the model. Under correct convergence,
    the maximized log-likelihood should therefore be (approximately) the
    same across all levels. A fit whose log-likelihood departs from the
    others by more than `tol_ll` has not found the global maximum.

2.  **Quantile monotonicity.** `mu_t^{(tau)}`, the conditional
    `tau`-quantile, must be non-decreasing in `tau` for every `t`, by
    definition of a quantile. This catches spurious optima the
    invariance check alone can miss (the likelihood surface can be
    nearly flat near a spurious optimum).

Fits flagged by either check are retried once from the parameter
estimates of each already-accepted level (not only the nearest one,
since the spurious basin can be reached from several "reasonable"
starting points), keeping the restart with the highest log-likelihood
among those that restore monotonicity with the levels already accepted.
This function does **not** cascade restarts further: an earlier
investigation found that a restart cascade can converge to a degenerate
region of the parameter space that, despite a nominally higher
log-likelihood, violates the invariance property itself (an unambiguous
sign of a spurious local optimum) – instead, unresolved cases are
flagged for manual inspection.

## See also

[`fit_copula()`](https://raydonal.github.io/kumarsurv/reference/fit_copula.md)
