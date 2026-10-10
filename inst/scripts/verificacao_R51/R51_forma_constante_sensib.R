## R51_forma_constante_sensib.R -- Rodada 51, item 51.14 (pedido de Pimentel, 2026-10-08).
## Rodar a partir de PACOTE_VALIDACAO/scripts:
##   LC_ALL=C.utf8 LANG=C.utf8 Rscript verificacao_R51/R51_forma_constante_sensib.R
## (1) Teste da razao de verossimilhancas dos tres coeficientes de forma sob AR <= 0.97,
##     com log-verossimilhancas nao arredondadas (S6.7 do suplementar).
## (2) Analise de sensibilidade: refaz a Tabela S12 (calib_controle_SEM.R) e a deteccao
##     pela CUSUM (02_farrington_SEM.R) com o modelo de forma constante (Kuma) e precisao
##     constante (beta), maximos restritos a AR <= 0.97 e, para comparacao, AR <= 0.996.
##     Os maximos restritos usam a mesma maquinaria de R50_tendencias_semestral.R.
## Os conjuntos de avaliacao (145 semanas sob controle, 18 do surto) sao os do artigo.
source("copula_ml.R")
dir_out <- file.path(DIR_RESULTADOS, "verificacao_R51"); dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
zz <- file(file.path(dir_out, "R51_forma_constante_sensib.txt"), open = "wt"); sink(zz, split = TRUE)
R <- readRDS("reais_result_sem.rds")
d <- R$d; X <- R$X; y <- R$y; phaseI <- R$phaseI; n <- nrow(d)
I <- seq_len(phaseI); II <- (phaseI + 1):n
Zc <- matrix(1, n, 1); kx <- ncol(X)
F_out <- readRDS("farrington_result.rds")$alarms_idx
C_out <- readRDS("cusum_result.rds")$flag
ctrl <- setdiff(II, union(F_out, C_out))
cat(sprintf("Fase II sob controle: %d semanas; surto (Farrington): %d semanas\n\n", length(ctrl), length(F_out)))

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
nll_I <- function(Zm, fam, tau) function(p) nll_copula(p, y[I], X[I, ], Zm[I, , drop = FALSE], 1, 1, fam, tau)

## ---------------------------------------------------------------- (1) razao de verossimilhancas
fk <- R$fits[["kuma 1 1"]]; k10 <- length(fk$par)
cat("=== (1) Teste RV dos coeficientes de forma (Kuma, tau = 0.5) ===\n")
cat(sprintf("modelo completo, ajuste do artigo: AR = %.4f, loglik = %.4f\n", fk$par[k10 - 1], fk$loglik))
sdk10 <- c(rep(1.5, kx), rep(.4, kx), .02, .05)
bfull <- modo(nll_I(R$Z, "kuma", 0.5), list(fk$par), k10 - 1, k10, sdk10, 0.97)
cat(sprintf("modelo completo, maximo com AR <= 0.97: AR = %.4f, loglik = %.4f\n", bfull$par[k10 - 1], -bfull$value))
f0 <- fit_copula(y[I], X[I, ], Zc[I, , drop = FALSE], p = 1, q = 1, family = "kuma", tau = 0.5); k7 <- length(f0$par)
st_int <- f0$par; st_int[k7 - 1] <- min(st_int[k7 - 1], 0.95)
sdk7 <- c(rep(1.5, kx), .4, .02, .05)
b0 <- modo(nll_I(Zc, "kuma", 0.5), list(st_int), k7 - 1, k7, sdk7, 0.97)
cat(sprintf("forma constante, maximo com AR <= 0.97: AR = %.4f, loglik = %.4f, AIC = %.2f\n",
            b0$par[k7 - 1], -b0$value, 2 * b0$value + 2 * k7))
LR <- 2 * (-bfull$value - (-b0$value))
cat(sprintf("RV = 2*(%.4f - %.4f) = %.4f, gl = %d, p = %.4f\n\n", -bfull$value, -b0$value, LR, k10 - k7,
            pchisq(LR, k10 - k7, lower.tail = FALSE)))

## ---------------------------------------------------------------- (2) sensibilidade
taus <- c(0.50, 0.75, 0.90, 0.95, 0.99)
pred <- function(par, family, tau = 0.5) {
  mu <- plogis(as.numeric(X %*% par[1:kx])); sh <- exp(par[kx + 1])
  ar <- par[kx + 2]; ma <- par[kx + 3]
  Ft <- if (family == "kuma") pkuma(y, mu, sh, tau) else pbetamp(y, mu, sh)
  eps <- qnorm(pmin(pmax(Ft, 1e-12), 1 - 1e-12))
  mom <- .one_step_moments(eps, .dl_recursion(as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = n - 1))))
  Q <- function(u) if (family == "kuma") qkuma(u, mu, sh, tau) else qbetamp(u, mu, sh)
  list(marg = function(t0) Q(rep(t0, n)), cond = function(t0) Q(pnorm(mom$m + sqrt(mom$s2) * qnorm(t0))),
       r = (eps - mom$m) / sqrt(mom$s2))
}
linha <- function(lim, rot, t0, ctrl_set = ctrl) data.frame(
  limite = rot, tau = t0,
  faseI = mean(y[I] <= lim[I]), faseII_ctrl = mean(y[ctrl_set] <= lim[ctrl_set]), faseII_todas = mean(y[II] <= lim[II]),
  surto_acima = sum(y[F_out] > lim[F_out]), lim_medio_I = mean(lim[I]))
cusum <- function(r, k = 0.5, h = 4) { cp <- numeric(n)
  for (i in II) cp[i] <- max(0, r[i] - k + cp[i - 1]); which(cp[II] > h) + phaseI }
descreve <- function(fl) if (length(fl) == 0) "nenhum alarme" else
  sprintf("%d semanas, de %s a %s; fora de 2023-W01..2023-W30: %d", length(fl), d$epiweek[min(fl)], d$epiweek[max(fl)],
          sum(!(d$epiweek[fl] >= "2023-W01" & d$epiweek[fl] <= "2023-W30")))
cat(sprintf("CUSUM do artigo (modelo completo, h = 4): %s\n", descreve(C_out)))

out <- list()
for (amax in c(0.97, 0.996)) {
  cat(sprintf("\n=== (2) Forma/precisao constante, maximos com AR <= %s ===\n", amax))
  fits <- list(); prev <- NULL
  for (t0 in taus) {
    fu <- try(fit_copula(y[I], X[I, ], Zc[I, , drop = FALSE], p = 1, q = 1, family = "kuma", tau = t0), silent = TRUE)
    st <- list(); if (!is.null(prev)) st <- c(st, list(prev))
    if (!inherits(fu, "try-error")) st <- c(st, list(fu$par))
    st <- c(st, list(b0$par))
    b <- modo(nll_I(Zc, "kuma", t0), st, k7 - 1, k7, sdk7, amax)
    fits[[as.character(t0)]] <- b; prev <- b$par
    cat(sprintf("  Kuma tau = %.2f: AR = %.4f  MA = %.4f  loglik = %.3f\n", t0, b$par[k7 - 1], b$par[k7], -b$value))
  }
  mu_mat <- sapply(fits, function(b) plogis(as.numeric(X %*% b$par[1:kx])))
  cat(sprintf("  quantis marginais nao monotonos em tau: %d de %d semanas\n",
              sum(apply(mu_mat, 1, function(r) any(diff(r) < -1e-8))), n))
  fb0 <- fit_copula(y[I], X[I, ], Zc[I, , drop = FALSE], p = 1, q = 1, family = "beta")
  stb <- fb0$par; stb[k7 - 1] <- min(stb[k7 - 1], amax - 0.01)
  bb <- modo(nll_I(Zc, "beta", 0.5), list(stb, fb0$par), k7 - 1, k7, sdk7, amax)
  cat(sprintf("  beta (precisao constante): AR = %.4f  loglik = %.3f  AIC = %.2f (irrestrito: AR = %.4f, AIC = %.2f)\n",
              bb$par[k7 - 1], -bb$value, 2 * bb$value + 2 * k7, fb0$par[k7 - 1], fb0$aic))
  cat(sprintf("  Kuma tau = 0.5: AIC = %.2f\n", 2 * fits[["0.5"]]$value + 2 * k7))
  pb <- pred(bb$par, "beta")
  pk5 <- pred(fits[["0.5"]]$par, "kuma", 0.5)
  fl <- cusum(pk5$r)
  cat(sprintf("  CUSUM (Kuma, forma constante, h = 4): %s\n", descreve(fl)))
  cat(sprintf("  semanas sinalizadas: %s\n", paste(d$epiweek[fl], collapse = ", ")))
  cat(sprintf("  em comum com a CUSUM do artigo: %d; so' nesta: %d; so' no artigo: %d\n",
              length(intersect(fl, C_out)), length(setdiff(fl, C_out)), length(setdiff(C_out, fl))))
  ctrl2 <- setdiff(II, union(F_out, fl))
  rows <- list()
  for (t0 in c(0.90, 0.95, 0.99)) {
    pk <- pred(fits[[as.character(t0)]]$par, "kuma", t0)
    rows <- c(rows, list(linha(pk$marg(t0), "Kuma marginal (ajuste em tau)", t0),
                         linha(pb$marg(t0), "Beta marginal (derivado da media)", t0),
                         linha(pk$cond(t0), "Kuma condicional (ajuste em tau)", t0),
                         linha(pb$cond(t0), "Beta condicional (derivado da media)", t0)))
  }
  tab <- do.call(rbind, rows)
  cat(sprintf("  Tabela S12 com forma/precisao constante (Fase II sob controle = as %d semanas do artigo):\n", length(ctrl)))
  options(width = 200); print(tab, digits = 3, row.names = FALSE)
  cat(sprintf("  Fase II sob controle redefinida com a CUSUM deste modelo: %d semanas; cobertura nelas:\n", length(ctrl2)))
  for (t0 in c(0.90, 0.95, 0.99)) {
    pk <- pred(fits[[as.character(t0)]]$par, "kuma", t0)
    cat(sprintf("    tau = %.2f: Kuma marginal %.3f, Kuma condicional %.3f, beta marginal %.3f, beta condicional %.3f\n", t0,
                mean(y[ctrl2] <= pk$marg(t0)[ctrl2]), mean(y[ctrl2] <= pk$cond(t0)[ctrl2]),
                mean(y[ctrl2] <= pb$marg(t0)[ctrl2]), mean(y[ctrl2] <= pb$cond(t0)[ctrl2])))
  }
  out[[as.character(amax)]] <- list(fits = fits, beta = bb, tab = tab, cusum = fl, ctrl2 = ctrl2)
}
sink(); close(zz)
saveRDS(list(LR = list(full = bfull, const = b0, LR = LR), sens = out, ctrl = ctrl, F_out = F_out, C_out = C_out),
        file.path(dir_out, "R51_forma_constante_sensib.rds"))
