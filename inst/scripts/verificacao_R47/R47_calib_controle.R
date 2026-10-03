## R47_calib_controle.R -- Rodada 47 (sandbox)
## Separa, na Fase II, a CALIBRACAO SOB CONTROLE (semanas fora do surto) da
## SENSIBILIDADE NO SURTO (semanas sinalizadas), para:
##   (a) limites MARGINAIS (os da tab:calib do v11): quantil marginal Kuma estimado
##       em tau (fit_copula_multi_tau) e quantil beta derivado do ajuste na media;
##   (b) limites CONDICIONAIS ao passado (quantil preditivo de um passo), obtidos
##       pela recursao de Durbin-Levinson com os parametros da Fase I fixos.
## Rotulos de surto: semanas sinalizadas pelo Farrington (farrington_result.rds)
## e pela CUSUM (cusum_result.rds).
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds")
d <- R$d; X <- R$X; Z <- R$Z; y <- R$y; phaseI <- R$phaseI; n <- nrow(d)
I <- seq_len(phaseI); II <- (phaseI + 1):n; yII <- y[II]
kx <- ncol(X); kz <- ncol(Z)
F_out <- readRDS("farrington_result.rds")$alarms_idx
C_out <- readRDS("cusum_result.rds")$flag
U_out <- sort(union(F_out, C_out))
ctrl  <- setdiff(II, U_out)
cat(sprintf("Fase II: %d semanas; Farrington %d; CUSUM %d; uniao %d; sob controle %d\n",
            length(II), length(F_out), length(C_out), length(U_out), length(ctrl)))
cat("semanas da uniao:", paste(d$epiweek[U_out], collapse = " "), "\n")

taus <- c(0.90, 0.95, 0.99)
fits <- fit_copula_multi_tau(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma",
                             taus = c(0.50, 0.75, 0.90, 0.95, 0.99), verbose = FALSE)
fk5 <- R$fits[["kuma 1 1"]]; fb <- R$fits[["beta 1 1"]]
cat(sprintf("tau do ajuste 'kuma 1 1' arquivado: %s\n", format(fk5$tau)))

## recursao preditiva com parametros fixos (Fase I) aplicada a toda a serie
pred <- function(par, family, tau = 0.5) {
  mu <- plogis(as.numeric(X %*% par[1:kx])); sh <- exp(as.numeric(Z %*% par[(kx + 1):(kx + kz)]))
  ar <- par[kx + kz + 1]; ma <- par[kx + kz + 2]
  Ft <- if (family == "kuma") pkuma(y, mu, sh, tau) else pbetamp(y, mu, sh)
  eps <- qnorm(pmin(pmax(Ft, 1e-12), 1 - 1e-12))
  rho <- as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = n - 1))
  mom <- .one_step_moments(eps, .dl_recursion(rho))
  Q <- function(u) if (family == "kuma") qkuma(u, mu, sh, tau) else qbetamp(u, mu, sh)
  list(mu = mu, sh = sh, m = mom$m, v = mom$s2,
       marg = function(t0) Q(rep(t0, n)),
       cond = function(t0) Q(pnorm(mom$m + sqrt(mom$s2) * qnorm(t0))))
}
pk5 <- pred(fk5$par, "kuma", 0.5); pb <- pred(fb$par, "beta")

rows <- list()
add <- function(lim, nome, tau0) {
  rows[[length(rows) + 1]] <<- data.frame(
    limite = nome, tau = tau0,
    calib_faseII = mean(y[II] <= lim[II]),
    calib_controle = mean(y[ctrl] <= lim[ctrl]),
    excede_controle = sum(y[ctrl] > lim[ctrl]),
    detecta_F = sum(y[F_out] > lim[F_out]), detecta_C = sum(y[C_out] > lim[C_out]))
}
for (t0 in taus) {
  ft <- fits[[as.character(t0)]]
  pkt <- pred(ft$par, "kuma", t0)
  add(pkt$marg(t0), "Kuma marginal, ajuste em tau (v11)", t0)
  add(pb$marg(t0),  "Beta marginal, derivado da media (v11)", t0)
  add(pk5$marg(t0), "Kuma marginal, derivado da mediana", t0)
  add(pkt$cond(t0), "Kuma condicional, ajuste em tau", t0)
  add(pk5$cond(t0), "Kuma condicional, derivado da mediana", t0)
  add(pb$cond(t0),  "Beta condicional, derivado da media", t0)
}
tab <- do.call(rbind, rows)
options(width = 160)
print(tab, digits = 3, row.names = FALSE)
cat(sprintf("\nalvo de excedencias sob controle (%d semanas): tau=0.90 -> %.1f; 0.95 -> %.1f; 0.99 -> %.1f\n",
            length(ctrl), 0.10 * length(ctrl), 0.05 * length(ctrl), 0.01 * length(ctrl)))
saveRDS(list(tab = tab, ctrl = ctrl, F_out = F_out, C_out = C_out), file.path(DIR_RESULTADOS, "R47_calib_controle.rds"))
