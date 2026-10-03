## perfil_ar1_exato.R -- Rodada 47.
## (1) Perfil de verossimilhanca do AR(1) da copula (Kumaraswamy e beta, Fase I), com
##     partidas nas duas bacias da verossimilhanca -- a do ajuste reportado e a vizinha
##     da fronteira de estacionariedade -- e grade ate 0.997. Usa a verossimilhanca de
##     copula_ml.R, que desde a Rodada 47 coincide com a exata junto a fronteira.
##     perfil_ar1.R (mantido) parte so' da bacia do ajuste reportado e para em 0.985.
## (2) Comparacao de periodicidade dos harmonicos com todos os ajustes sob a mesma
##     restricao: o maximo com AR <= 0.97 e o maximo com AR <= 0.996.
## Maximos restritos: AR = -0.997 + (amax + 0.997) plogis(a) e MA = tanh(b), otimizados
## sem restricao em (a, b). Uma barreira (penalidade 1e10 fora da faixa) faz o BFGS
## parar antes do maximo restrito, e o resultado passa a depender do caminho do otimizador.
source("copula_ml.R")
out <- file.path(DIR_RESULTADOS, "resultado_perfil_ar1_exato.txt"); zz <- file(out, "wt"); sink(zz, split = TRUE)
R <- readRDS("reais_result_sem.rds"); d <- R$d; phaseI <- R$phaseI; I <- seq_len(phaseI)
y <- R$y[I]; X <- R$X[I, ]; Z <- R$Z[I, ]
## maximo de -nll com AR (posicao iar) em [-0.997, amax] e MA (posicao ima) em (-1, 1)
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
cat("=== (1) PERFIL DO AR(1) DA COPULA, FASE I ===\n")
sdv <- c(1.5, 4, .5, .5, .4, .8, .1, .1)
for (fam in c("kuma", "beta")) {
  f <- R$fits[[paste(fam, "1 1")]]; k <- length(f$par)
  nll <- function(p) nll_copula(p, y, X, Z, 1, 1, fam, 0.5)
  obj <- function(p) if (abs(p[10]) >= 1) 1e10 else nll(p)
  p_front <- f$par; p_front[9] <- 0.99
  b <- modo(nll, list(f$par, p_front), 9, 10, c(sdv, .02, .05), 0.997)
  cat(sprintf("\n[%s] ajuste reportado: loglik %.3f, AR %.4f | maximo com AR <= 0.997: %.3f, AR %.4f\n",
              fam, f$loglik, f$par[9], -b$value, b$par[9]))
  pf <- function(a) { bl <- -Inf; for (s in list(f$par[-9], b$par[-9])) {
      nf <- function(q) { p <- numeric(k); p[-9] <- q; p[9] <- a; obj(p) }
      o <- try(optim(s, nf, method = "BFGS", control = list(maxit = 3000, reltol = 1e-11)), silent = TRUE)
      if (!inherits(o, "try-error") && -o$value > bl) bl <- -o$value }; bl }
  grid <- c(seq(0.85, 0.98, by = 0.01), 0.985, 0.99, 0.993, 0.995, 0.997)
  lp <- sapply(grid, pf); lmax <- max(lp)
  print(round(cbind(AR = grid, loglik = lp, dif_max = lmax - lp, dif_reportado = f$loglik - lp), 3))
  cat(sprintf("[%s] maximo do perfil em AR = %.3f; o perfil %s na fronteira da grade\n", fam, grid[which.max(lp)],
              if (which.max(lp) == length(grid)) "cresce ate'" else "nao cresce ate'"))
}
cat("\n=== (2) PERIODICIDADE DOS HARMONICOS SOB A MESMA RESTRICAO (Kuma, forma so' com intercepto) ===\n")
y0 <- ifelse(d$incidence[I] == 0, 1e-7, d$incidence[I]); s <- I; tr <- (s - mean(s)) / 100
per <- c("52 (anual)" = 52, "26 (semestral)" = 26, "13 (trimestral)" = 13, "4.33 (mensal)" = 365.25 / 12 / 7, "sem harmonicos" = NA)
tab <- data.frame()
for (nm in names(per)) {
  P <- per[[nm]]; Xp <- if (is.na(P)) cbind(1, tr) else cbind(1, tr, cos(2 * pi * s / P), sin(2 * pi * s / P)); Zp <- matrix(1, phaseI, 1)
  f0 <- fit_copula(y0, Xp, Zp, p = 1, q = 1, family = "kuma", tau = 0.5); k <- length(f0$par)
  nll <- function(p) nll_copula(p, y0, Xp, Zp, 1, 1, "kuma", 0.5)
  sdk <- c(rep(1.5, ncol(Xp)), .4, .02, .05)
  st_int <- f0$par; st_int[k - 1] <- min(st_int[k - 1], 0.95); st_fr <- f0$par; st_fr[k - 1] <- 0.99
  bi <- modo(nll, list(st_int), k - 1, k, sdk, 0.97); bf <- modo(nll, list(st_fr, st_int), k - 1, k, sdk, 0.996)
  tab <- rbind(tab, data.frame(periodo = nm, npar = k, fit_copula_AR = round(f0$par[k - 1], 3), fit_copula_AIC = round(-2 * f0$loglik + 2 * k, 2),
                               AR_ate_0.97 = round(bi$par[k - 1], 3), AIC_ate_0.97 = round(2 * bi$value + 2 * k, 2),
                               AR_ate_0.996 = round(bf$par[k - 1], 3), AIC_ate_0.996 = round(2 * bf$value + 2 * k, 2)))
}
options(width = 200); print(tab, row.names = FALSE)
cat("(AR igual a 0.970 ou 0.996 indica maximo na borda da restricao)\n")
sink(); close(zz)
saveRDS(tab, file.path(DIR_RESULTADOS, "perfil_ar1_exato_periodicidade.rds"))
