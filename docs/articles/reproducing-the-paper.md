# Reproducing the paper

This package exports the model-fitting functions behind *Copula-Based
Kumaraswamy Regression for Anomaly Detection in Doubly Bounded Time
Series* (Ospina, Pimentel and Cribari-Neto). Every table, figure and
numerical value in the paper and its supplementary material is generated
by a script archived in `inst/scripts/` (see the README’s
result-to-script table); this vignette walks through the core
model-fitting workflow directly with the exported functions, on the
paper’s own data.

## The data

``` r
library(kumarsurv)
data(platina)
head(platina)
#>    epiweek cases   incidence       sin       cos   trend
#> 1 2017-W01     0 0.000000000 0.1205367 0.9927089 -0.2345
#> 2 2017-W02     0 0.000000000 0.2393157 0.9709418 -0.2335
#> 3 2017-W03     0 0.000000000 0.3546049 0.9350162 -0.2325
#> 4 2017-W04     1 0.003305785 0.4647232 0.8854560 -0.2315
#> 5 2017-W05     2 0.006611570 0.5680647 0.8229839 -0.2305
#> 6 2017-W06     1 0.003305785 0.6631227 0.7485107 -0.2295
```

`platina$incidence` is `10 * cases / population`: a scaling that only
places the series in `(0, 1)`, unrelated to the per-ten-thousand rate
used to select the municipality (Section 6.1 of the paper). `48.9%` of
weeks have no recorded visits; since the continuous marginal
distributions considered here have open support in `(0, 1)`, zeros are
replaced by a boundary value for estimation only:

``` r
y <- platina$incidence
y[y == 0] <- 1e-7
```

## Fitting the model

Location and shape submodels with a trend and the first annual harmonic,
Gaussian copula with ARMA(1,1) dependence, Kumaraswamy margin at the
median (`tau = 0.5`):

``` r
X <- cbind(intercept = 1, trend = platina$trend, sin = platina$sin, cos = platina$cos)
Z <- cbind(intercept = 1, trend = platina$trend, sin = platina$sin, cos = platina$cos)
fit <- fit_copula(y, X, Z, p = 1, q = 1, family = "kuma", tau = 0.5)
fit$par
#>    mu.intercept        mu.trend          mu.sin          mu.cos shape.intercept 
#>     -9.17613580      7.65538429      0.05949805     -0.12771616     -1.64090693 
#>     shape.trend       shape.sin       shape.cos             ar1             ma1 
#>      0.65964750     -0.03409768     -0.01226825      0.91389214     -0.71817252
fit$loglik
#> [1] 3516.435
```

## Predictive quantile residuals

``` r
r <- qresiduals(fit)$r
summary(r)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -2.061949 -0.853502 -0.039043 -0.004756  0.883398  2.528820
```

Under correct specification, these are approximately i.i.d. standard
normal (the paper’s Proposition 1); they feed the goodness-of-fit tests,
the one-step prediction bands and the CUSUM chart
(`inst/scripts/13_bootstrap_gof.R`, `09_efeito_correlacao.R`,
`11_arl.R`, `02_farrington_SEM.R`).

## Fitting at several quantile levels

**Do not call
[`fit_copula()`](https://raydonal.github.io/kumarsurv/reference/fit_copula.md)
directly at `tau != 0.5`** without also checking the diagnostics
[`fit_copula_multi_tau()`](https://raydonal.github.io/kumarsurv/reference/fit_copula_multi_tau.md)
runs: the likelihood surface can have more than one interior stationary
point at a given `tau` (the paper’s Section 6.3/Table 3 and pitfall 1),
and a single-start fit is not guaranteed to land on the one with the
higher likelihood.

``` r
fits <- fit_copula_multi_tau(y, X, Z, p = 1, q = 1, family = "kuma",
                              taus = c(0.50, 0.90, 0.95), verbose = FALSE)
attr(fits, "diagnostico")
#>       tau   loglik  fonte monotono invariante_ll aprovado
#> 0.5  0.50 3516.435 padrao     TRUE          TRUE     TRUE
#> 0.9  0.90 3516.503 padrao     TRUE          TRUE     TRUE
#> 0.95 0.95 3516.514 padrao     TRUE          TRUE     TRUE
```

The `aprovado` column flags levels that passed both the
quantile-monotonicity and log-likelihood-invariance checks; a `FALSE`
means the fit needs manual inspection before its quantile is used as a
control limit.

## Simulating from the fitted model

Used for Monte Carlo validation and for the parametric bootstrap behind
the paper’s goodness-of-fit tests and prediction bands:

``` r
set.seed(2026)
y_sim <- sim_copula(X, Z, beta = fit$par[1:4], gamma = fit$par[5:8],
                     ar = fit$par["ar1"], ma = fit$par["ma1"],
                     family = "kuma", tau = 0.5)
head(y_sim)
#> Time Series:
#> Start = 1 
#> End = 6 
#> Frequency = 1 
#> [1] 9.275259e-08 2.859146e-06 2.623463e-07 7.664494e-04 4.918441e-09
#> [6] 2.164632e-05
```
