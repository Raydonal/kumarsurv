# Weekly acute diarrheal disease incidence in Platina, SP (semiannual harmonic)

Same series as
[platina](https://raydonal.github.io/kumarsurv/reference/platina.md),
with an additional first semiannual harmonic (period 26 weeks), the
periodic component actually adopted in the paper (Section 6.2: a
semiannual period gives the smallest AIC, consistent with two annual
transmission peaks being epidemiologically plausible for ADD).

## Usage

``` r
platina_sem
```

## Format

A data frame with 470 rows and 8 columns:

- epiweek:

  Character, ISO epidemiological week, e.g. `"2017-W01"`.

- cases:

  Integer, number of visits recorded in the week.

- incidence:

  Numeric, `10 * cases / population`, in `(0, 1)`.

- sin, cos:

  Numeric, first annual harmonic (period 52 weeks).

- trend:

  Numeric, linear time trend, centered and scaled.

- sin26, cos26:

  Numeric, first semiannual harmonic (period 26 weeks).

## Source

Primary Health Care records, Brazilian Ministry of Health, 2017-2025.
