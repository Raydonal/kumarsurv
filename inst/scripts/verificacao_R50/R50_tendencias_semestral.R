## R50_tendencias_semestral.R -- Rodada 50 (verificacao; rodar por rodar_verificacao_R50.R).
## A Parte 2 de analise_covariaveis.R compara especificacoes de tendencia com os
## harmonicos ANUAIS (d$cos, d$sin), mas o modelo adotado e' semestral. Refaz a
## comparacao com os harmonicos semestrais, forma so' com intercepto, (i) por
## fit_copula irrestrito, como a Parte 2, e (ii) pelos maximos restritos a
## AR <= 0.97 e AR <= 0.996, com a mesma maquinaria de perfil_ar1_exato.R.
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); d <- R$d; phaseI <- R$phaseI; I <- seq_len(phaseI)
y0 <- ifelse(d$incidence[I] == 0, 1e-7, d$incidence[I]); s <- I
modo <- function(nll, st_list, iar, ima, sdk, amax, n0 = 12, seed = 47) {
  lo <- -0.997
  to_p <- function(u) { p <- u; p[iar] <- lo + (amax - lo) * plogis(u[iar]); p[ima] <- tanh(u[ima]); p }
  to_u <- function(p) { u <- p; a <- min(max(p[iar], lo + 1e-6), amax - 1e-6)
    u[iar] <- qlogis((a - lo) / (amax - lo)); u[ima] <- atanh(min(max(p[ima], -0.999), 0.999)); u }
  set.seed(seed); starts <- c(st_list, lapply(seq_len(n0), function(i) st_list[[1]] + rnorm(length(st_list[[1]]), 0, sdk)))
  best <- NULL
  for (p0 in starts) {
    o <- try(optim(to_u(p0), function(u) nll(to_p(u)), method = "BFGS", control = list(maxit = 5000, reltol = 1e-12)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- list(value = o$value, par = to_p(o$par))
  }
  best
}
cs <- d$cos26[I]; sn <- d$sin26[I]
trends <- list("linear s" = (s - mean(s)) / 100,
               "sqrt(s)"  = { v <- sqrt(s); (v - mean(v)) / sd(v) },
               "log(s)"   = { v <- log(s); (v - mean(v)) / sd(v) },
               "sem tend." = NULL)
tab <- data.frame()
for (nm in names(trends)) {
  tr <- trends[[nm]]
  Xp <- if (is.null(tr)) cbind(1, cs, sn) else cbind(1, tr, cs, sn); Zp <- matrix(1, phaseI, 1)
  f0 <- fit_copula(y0, Xp, Zp, p = 1, q = 1, family = "kuma", tau = 0.5); k <- length(f0$par)
  nll <- function(p) nll_copula(p, y0, Xp, Zp, 1, 1, "kuma", 0.5)
  sdk <- c(rep(1.5, ncol(Xp)), .4, .02, .05)
  st_int <- f0$par; st_int[k - 1] <- min(st_int[k - 1], 0.95); st_fr <- f0$par; st_fr[k - 1] <- 0.99
  bi <- modo(nll, list(st_int), k - 1, k, sdk, 0.97); bf <- modo(nll, list(st_fr, st_int), k - 1, k, sdk, 0.996)
  tab <- rbind(tab, data.frame(tendencia = nm, npar = k,
    livre_AR = round(f0$par[k - 1], 3), livre_ll = round(f0$loglik, 2), livre_AIC = round(-2 * f0$loglik + 2 * k, 2),
    AR_097 = round(bi$par[k - 1], 3), ll_097 = round(-bi$value, 2), AIC_097 = round(2 * bi$value + 2 * k, 2),
    AR_0996 = round(bf$par[k - 1], 3), ll_0996 = round(-bf$value, 2), AIC_0996 = round(2 * bf$value + 2 * k, 2)))
  cat(nm, "ok\n")
}
options(width = 220); print(tab, row.names = FALSE)
saveRDS(tab, file.path(DIR_RESULTADOS, "verificacao_R50", "R50_tendencias_semestral.rds"))
