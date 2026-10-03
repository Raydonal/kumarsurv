#' Weekly acute diarrheal disease incidence in Platina, SP (annual harmonic)
#'
#' Weekly records of acute diarrheal disease (ADD) in the municipality of
#' Platina, Sao Paulo state, Brazil (population 3025), epidemiological week
#' 2017-W01 to 2025-W53 (`n = 470`). The first 300 observations are Phase I
#' (estimation and calibration) in the paper; the remaining 170 are Phase II
#' (monitoring). `incidence` is `10 * cases / population`, a scaling that only
#' places the series in `(0, 1)`, unrelated to the per-ten-thousand rate used
#' to select the municipality. `sin`/`cos` are the first annual harmonic.
#'
#' @format A data frame with 470 rows and 6 columns:
#' \describe{
#'   \item{epiweek}{Character, ISO epidemiological week, e.g. `"2017-W01"`.}
#'   \item{cases}{Integer, number of visits recorded in the week.}
#'   \item{incidence}{Numeric, `10 * cases / population`, in `(0, 1)`.}
#'   \item{sin, cos}{Numeric, first annual harmonic (period 52 weeks).}
#'   \item{trend}{Numeric, linear time trend, centered and scaled.}
#' }
#' @source Primary Health Care records, Brazilian Ministry of Health, 2017-2025.
"platina"

#' Weekly acute diarrheal disease incidence in Platina, SP (semiannual harmonic)
#'
#' Same series as [platina], with an additional first semiannual harmonic
#' (period 26 weeks), the periodic component actually adopted in the paper
#' (Section 6.2: a semiannual period gives the smallest AIC, consistent with
#' two annual transmission peaks being epidemiologically plausible for ADD).
#'
#' @format A data frame with 470 rows and 8 columns:
#' \describe{
#'   \item{epiweek}{Character, ISO epidemiological week, e.g. `"2017-W01"`.}
#'   \item{cases}{Integer, number of visits recorded in the week.}
#'   \item{incidence}{Numeric, `10 * cases / population`, in `(0, 1)`.}
#'   \item{sin, cos}{Numeric, first annual harmonic (period 52 weeks).}
#'   \item{trend}{Numeric, linear time trend, centered and scaled.}
#'   \item{sin26, cos26}{Numeric, first semiannual harmonic (period 26 weeks).}
#' }
#' @source Primary Health Care records, Brazilian Ministry of Health, 2017-2025.
"platina_sem"
