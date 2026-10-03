# Kumaraswamy distribution function, quantile-parametrized

Kumaraswamy distribution function, quantile-parametrized

## Usage

``` r
pkuma(y, mu, theta, tau = 0.5)
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

## Value

Numeric vector, `F(y)`.

## See also

[`dkuma()`](https://raydonal.github.io/kumarsurv/reference/dkuma.md),
[`qkuma()`](https://raydonal.github.io/kumarsurv/reference/qkuma.md)
