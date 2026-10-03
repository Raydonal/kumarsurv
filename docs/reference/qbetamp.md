# Beta quantile function, mean-precision parametrization

Unlike
[`qkuma()`](https://raydonal.github.io/kumarsurv/reference/qkuma.md),
this requires the numerical inversion built into
[`stats::qbeta()`](https://rdrr.io/r/stats/Beta.html): the beta
distribution function has no closed form.

## Usage

``` r
qbetamp(u, mu, phi)
```

## Arguments

- u:

  Numeric vector of probabilities in `(0, 1)`.

- mu:

  Numeric vector, the mean, in `(0, 1)`.

- phi:

  Numeric vector, the precision parameter, `> 0`.

## Value

Numeric vector, `F^{-1}(u)`.
