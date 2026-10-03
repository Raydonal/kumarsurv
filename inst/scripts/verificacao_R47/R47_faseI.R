source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); d <- R$d; X <- R$X; Z <- R$Z; y <- R$y; phaseI <- R$phaseI; n <- nrow(d)
I <- seq_len(phaseI); II <- (phaseI + 1):n; kx <- ncol(X)
F_out <- readRDS("farrington_result.rds")$alarms_idx; C_out <- readRDS("cusum_result.rds")$flag
ctrl <- setdiff(II, union(F_out, C_out)); yr <- d$incidence
cat("quantis empiricos da incidencia (Fase I):", sprintf("%.4f", quantile(yr[I], c(.5,.75,.9,.95,.99), type = 7)), "\n")
cat("quantis empiricos (Fase II sob controle):", sprintf("%.4f", quantile(yr[ctrl], c(.5,.75,.9,.95,.99), type = 7)), "\n")
cat("1 caso = ", sprintf("%.5f", 10/3025), "; max Fase I =", max(yr[I]), " max Fase II controle =", max(yr[ctrl]), "\n")
fits <- fit_copula_multi_tau(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma", taus = c(0.50,0.75,0.90,0.95,0.99), verbose = FALSE)
fb <- R$fits[["beta 1 1"]]; mub <- plogis(as.numeric(X %*% fb$par[1:kx])); phib <- exp(as.numeric(Z %*% fb$par[(kx+1):(2*kx)]))
for (t0 in c(0.90, 0.95, 0.99)) { muk <- plogis(as.numeric(X %*% fits[[as.character(t0)]]$par[1:kx])); qb <- qbeta(t0, mub*phib, (1-mub)*phib)
  cat(sprintf("tau=%.2f  cobertura Fase I (dentro da amostra): Kuma %.3f  Beta %.3f | em casos: limite medio Fase I Kuma %.1f  Beta %.1f\n",
      t0, mean(y[I] <= muk[I]), mean(y[I] <= qb[I]), mean(muk[I])*302.5, mean(qb[I])*302.5)) }
