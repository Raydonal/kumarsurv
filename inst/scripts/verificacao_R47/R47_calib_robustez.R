## R47_calib_robustez.R -- Rodada 47 (sandbox): robustez da avaliacao sob controle
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); d <- R$d; y <- R$y; phaseI <- R$phaseI; n <- nrow(d)
I <- seq_len(phaseI); II <- (phaseI + 1):n
F_out <- readRDS("farrington_result.rds")$alarms_idx; C_out <- readRDS("cusum_result.rds")$flag
sets <- list("Fase II inteira" = II, "sem Farrington" = setdiff(II, F_out),
             "sem Farrington+CUSUM" = setdiff(II, union(F_out, C_out)))
cat("Fase II comeca em", d$epiweek[phaseI + 1], "e termina em", d$epiweek[n], "\n")
especs <- list("com tendencia (v11)" = cbind(1, d$trend, d$cos26, d$sin26),
               "sem tendencia"       = cbind(1, d$cos26, d$sin26))
for (nm in names(especs)) {
  X <- especs[[nm]]; Z <- X; kx <- ncol(X)
  fits <- fit_copula_multi_tau(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma",
                               taus = c(0.50, 0.75, 0.90, 0.95, 0.99), verbose = FALSE)
  fb <- fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "beta")
  mub <- plogis(as.numeric(X %*% fb$par[1:kx])); phib <- exp(as.numeric(Z %*% fb$par[(kx + 1):(2 * kx)]))
  cat(sprintf("\n=== %s: loglik Kuma(0.5)=%.2f  beta=%.2f  AIC Kuma=%.2f beta=%.2f ===\n", nm,
              fits[["0.5"]]$loglik, if (!is.null(fb$loglik)) fb$loglik else NA, NA, NA))
  for (t0 in c(0.90, 0.95, 0.99)) {
    muk <- plogis(as.numeric(X %*% fits[[as.character(t0)]]$par[1:kx])); qb <- qbeta(t0, mub * phib, (1 - mub) * phib)
    s <- sapply(sets, function(ix) sprintf("K %.3f / B %.3f", mean(y[ix] <= muk[ix]), mean(y[ix] <= qb[ix])))
    cat(sprintf("tau=%.2f  %s\n", t0, paste(sprintf("%s: %s", names(s), s), collapse = " | ")))
    cat(sprintf("          limite medio Fase I: K %.4f B %.4f | Fase II: K %.4f B %.4f | excedencias surto(F): K %d B %d\n",
                mean(muk[I]), mean(qb[I]), mean(muk[II]), mean(qb[II]), sum(y[F_out] > muk[F_out]), sum(y[F_out] > qb[F_out])))
  }
}
