###############################################################################
## sensib_zeros_controle.R -- Rodada 47. Sensibilidade a epsilon (valor que
## substitui os zeros), com a calibracao do limite de tau = 0.95 medida SOB
## CONTROLE (Fase I e Fase II fora do surto) e a deteccao nas 18 semanas do
## surto sinalizadas pelo Farrington, como em calib_controle_SEM.R. Os ajustes
## sao os mesmos de analise_sensibilidade_zeros.R (que mede a calibracao em
## todas as 170 semanas da Fase II e continua registrado no pacote).
###############################################################################
source("copula_ml.R")
d <- read.csv("dda_platina_sem.csv", stringsAsFactors = FALSE)
y_raw <- d$incidence
X <- as.matrix(cbind(1, d$trend, d$cos26, d$sin26)); Z <- X
phaseI <- 300; I <- seq_len(phaseI); n <- nrow(d); II <- (phaseI + 1):n
F_out <- readRDS("farrington_result.rds")$alarms_idx; C_out <- readRDS("cusum_result.rds")$flag
ctrl <- setdiff(II, union(F_out, C_out))
out <- file.path(DIR_RESULTADOS, "resultado_sensib_zeros_controle.txt"); zz <- file(out, "wt"); sink(zz, split = TRUE)
cat("=== SENSIBILIDADE A epsilon: calibracao do limite de tau = 0.95 sob controle ===\n")
cat(sprintf("Fase I: %d semanas; Fase II sob controle: %d semanas; surto (Farrington): %d semanas\n\n", length(I), length(ctrl), length(F_out)))
cat(sprintf("%-7s | %-9s %-9s | %-21s | %-21s | %-21s\n", "eps", "AIC Kuma", "AIC beta", "Fase I (K | B)", "Fase II ctrl (K | B)", "surto acima (K | B)"))
res <- list()
for (eps in c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)) {
  y <- ifelse(y_raw == 0, eps, y_raw)
  fk <- fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma", tau = 0.5)
  fb <- fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "beta", tau = 0.5)
  fits <- fit_copula_multi_tau(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma",
                               taus = c(0.50, 0.75, 0.90, 0.95, 0.99), verbose = FALSE)
  muk <- plogis(as.numeric(X %*% fits[["0.95"]]$par[1:ncol(X)]))
  mub <- plogis(as.numeric(X %*% fb$par[1:ncol(X)])); phib <- exp(as.numeric(Z %*% fb$par[(ncol(X) + 1):(2 * ncol(X))]))
  qb <- qbeta(0.95, mub * phib, (1 - mub) * phib)
  r <- c(eps = eps, kaic = fk$aic, baic = fb$aic,
         kI = mean(y_raw[I] <= muk[I]), bI = mean(y_raw[I] <= qb[I]),
         kC = mean(y_raw[ctrl] <= muk[ctrl]), bC = mean(y_raw[ctrl] <= qb[ctrl]),
         kII = mean(y_raw[II] <= muk[II]), bII = mean(y_raw[II] <= qb[II]),
         kS = sum(y_raw[F_out] > muk[F_out]), bS = sum(y_raw[F_out] > qb[F_out]))
  res[[length(res) + 1]] <- r
  cat(sprintf("%-7.0e | %9.2f %9.2f | %.3f | %.3f         | %.3f | %.3f         | %2d | %2d      (Fase II inteira: %.3f | %.3f)\n",
              eps, fk$aic, fb$aic, r["kI"], r["bI"], r["kC"], r["bC"], r["kS"], r["bS"], r["kII"], r["bII"]))
}
sink(); close(zz)
saveRDS(do.call(rbind, res), file.path(DIR_RESULTADOS, "sensib_zeros_controle.rds"))
