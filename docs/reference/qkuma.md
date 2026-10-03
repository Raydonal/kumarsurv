# Kumaraswamy quantile function, quantile-parametrized

Closed form: no numerical inversion is needed, unlike the beta
distribution. This is the property the paper's copula construction,
predictive quantile residuals and prediction bands rely on throughout.

## Usage

``` r
qkuma(u, mu, theta, tau = 0.5)
```

## Arguments

- u:

  Numeric vector of probabilities in `(0, 1)`.

- mu:

  Numeric vector, the `tau`-quantile, in `(0, 1)`.

- theta:

  Numeric vector, the shape parameter, `> 0`.

- tau:

  Quantile level, in `(0, 1)`. Default `0.5` (median).

## Value

Numeric vector, `F^{-1}(u)`.

## See also

[`dkuma()`](https://raydonal.github.io/kumarsurv/reference/dkuma.md),
[`pkuma()`](https://raydonal.github.io/kumarsurv/reference/pkuma.md)
