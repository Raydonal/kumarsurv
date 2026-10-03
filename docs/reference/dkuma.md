# Kumaraswamy density, quantile-parametrized

Density of the Kumaraswamy distribution reparametrized so that `mu` is
its `tau`-quantile, `mu = F^{-1}(tau)`, rather than one of the two
standard shape parameters. For `tau = 0.5`, `mu` is the median, the
parametrization of Mitnik and Baek (2013).

## Usage

``` r
dkuma(y, mu, theta, tau = 0.5, log = FALSE)
```

## Arguments

- y:

  Numeric vector in `(0, 1)`.

- mu:

  Numeric vector, the `tau`-quantile, in `(0, 1)`.

- theta:

  Numeric vector, the shape parameter, `> 0`.

- tau:

  Quantile level, in `(0, 1)`. Default `0.5` (median).

- log:

  Logical; if `TRUE`, the log-density is returned.

## Value

Numeric vector of (log-)densities.

## See also

[`pkuma()`](https://raydonal.github.io/kumarsurv/reference/pkuma.md),
[`qkuma()`](https://raydonal.github.io/kumarsurv/reference/qkuma.md)
