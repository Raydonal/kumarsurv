## R47_periodo_modos.R -- comparacao de periodicidade: em que modo (AR) caiu cada ajuste de
## fit_copula (analise_harmonicos_inflacao.R), e partidas multiplas com |AR| <= 0.996.
## Diagnostico da auditoria. A coluna "melhor" usa uma penalidade fora da faixa, que faz o
## BFGS parar antes do maximo restrito: ela SUBESTIMA o maximo e depende do caminho do
## otimizador. Os maximos restritos da Tabela S12 vem de perfil_ar1_exato.R, que reparametriza
## o AR e o MA e nao tem esse problema.
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); d <- R$d; phaseI <- R$phaseI
y <- ifelse(d$incidence[1:phaseI] == 0, 1e-7, d$incidence[1:phaseI]); s <- seq_len(phaseI); tr <- (s - mean(s)) / 100
per <- c("52 (anual)" = 52, "26 (semestral)" = 26, "13 (trimestral)" = 13, "4.33 (mensal)" = 365.25 / 12 / 7, "sem harmonicos" = NA)
for (nm in names(per)) {
  P <- per[[nm]]; X <- if (is.na(P)) cbind(1, tr) else cbind(1, tr, cos(2 * pi * s / P), sin(2 * pi * s / P)); Z <- matrix(1, phaseI, 1)
  f <- fit_copula(y, X = X, Z = Z, family = "kuma", p = 1, q = 1, tau = 0.5); k <- length(f$par)
  obj <- function(p) { if (abs(p[k - 1]) > 0.996 || abs(p[k]) >= 1) return(1e10); nll_copula(p, y, X, Z, 1, 1, "kuma", 0.5) }
  set.seed(26); st <- c(list(f$par), lapply(1:15, function(i) f$par + rnorm(k, 0, c(rep(1.5, ncol(X)), 0.4, 0.03, 0.06))))
  st <- c(st, lapply(c(0.90, 0.95, 0.99), function(a) { p <- f$par; p[k - 1] <- a; p }))
  o <- lapply(st, function(p0) try(optim(p0, obj, method = "BFGS", control = list(maxit = 3000, reltol = 1e-11)), silent = TRUE))
  o <- o[!sapply(o, inherits, "try-error")]; ll <- sapply(o, function(z) -z$value); b <- o[[which.max(ll)]]
  modos <- unique(round(cbind(ll, sapply(o, function(z) z$par[k - 1])), 2)); modos <- modos[order(-modos[, 1]), , drop = FALSE]
  cat(sprintf("%-16s fit_copula: loglik %.2f AR %.3f AIC %.2f | melhor (|AR|<=0.996): loglik %.2f AR %.3f AIC %.2f | modos: %s\n",
              nm, f$loglik, f$par[k - 1], -2 * f$loglik + 2 * k, -b$value, b$par[k - 1], 2 * b$value + 2 * k,
              paste(sprintf("%.2f@%.3f", modos[1:min(4, nrow(modos)), 1], modos[1:min(4, nrow(modos)), 2]), collapse = " ")))
}
