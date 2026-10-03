## Code to prepare the `platina` dataset, kept per usethis convention
## (data-raw/, not installed with the package; run once to regenerate data/platina.rda).

platina <- read.csv("data-raw/dda_platina.csv", stringsAsFactors = FALSE)
platina_sem <- read.csv("data-raw/dda_platina_sem.csv", stringsAsFactors = FALSE)

usethis::use_data(platina, overwrite = TRUE)
usethis::use_data(platina_sem, overwrite = TRUE)
