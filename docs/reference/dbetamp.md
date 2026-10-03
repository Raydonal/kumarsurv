# Beta density, mean-precision parametrization

The mean-precision beta regression parametrization of Ferrari and
Cribari-Neto (2004), used throughout as the comparison marginal
distribution.

## Usage

``` r
dbetamp(y, mu, phi, log = FALSE)
```

## Arguments

- y:

  Numeric vector in `(0, 1)`.

- mu:

  Numeric vector, the mean, in `(0, 1)`.

- phi:

  Numeric vector, the precision parameter, `> 0`.

- log:

  Logical; if `TRUE`, the log-density is returned.

## Value

Numeric vector of (log-)densities.

## See also

[`pbetamp()`](https://raydonal.github.io/kumarsurv/reference/pbetamp.md),
[`qbetamp()`](https://raydonal.github.io/kumarsurv/reference/qbetamp.md)
