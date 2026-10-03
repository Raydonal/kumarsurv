# Weekly acute diarrheal disease incidence in Platina, SP (annual harmonic)

Weekly records of acute diarrheal disease (ADD) in the municipality of
Platina, Sao Paulo state, Brazil (population 3025), epidemiological week
2017-W01 to 2025-W53 (`n = 470`). The first 300 observations are Phase I
(estimation and calibration) in the paper; the remaining 170 are Phase
II (monitoring). `incidence` is `10 * cases / population`, a scaling
that only places the series in `(0, 1)`, unrelated to the
per-ten-thousand rate used to select the municipality. `sin`/`cos` are
the first annual harmonic.

## Usage

``` r
platina
```

## Format

A data frame with 470 rows and 6 columns:

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

## Source

Primary Health Care records, Brazilian Ministry of Health, 2017-2025.
