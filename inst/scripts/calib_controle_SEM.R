## calib_controle_SEM.R -- Rodada 47: calibracao dos limites de controle AVALIADA
## SOB CONTROLE, e deteccao nas semanas do surto -- harmonico SEMESTRAL.
##
## Motivo: calib_SEM.R (tab:calib ate o main_v11) mede a fracao de observacoes da
## Fase II abaixo do limite em TODAS as 170 semanas, inclusive as do surto de 2023.
## Semanas de surto excedem o limite por construcao do problema, de modo que essa
## fracao mistura alarmes falsos com alarmes verdadeiros: um limite alto demais sob
## controle e que deixa de sinalizar parte do surto pode parecer "calibrado". A
## calibracao de um limite de controle se afere nas semanas sob controle; o surto
## afere a sensibilidade.
##
## Conjuntos de avaliacao:
##   Fase I      : as n_I = 300 semanas de ajuste (nenhuma anomalia sinalizada nelas);
##   Fase II ctrl: semanas da Fase II fora do surto, isto e, nao sinalizadas pelo
##                 Farrington nem pela CUSUM (170 - 25 = 145 semanas);
##   surto       : as 18 semanas sinalizadas pelo Farrington (referencia externa).
## Limites:
##   marginais    : Kuma com a distribuicao marginal ajustada diretamente em tau
##                  (fit_copula_multi_tau, como em calib_SEM.R) e beta derivada do
##                  ajuste na media (os limites da tab:calib do main_v11);
##   condicionais : tau-quantil preditivo de um passo, F_t^{-1}(Phi(m + sqrt(v) z_tau)),
##                  com os parametros da Fase I fixos (quantil condicional ao passado,
##                  que e' o U_t da eq. (limiar) do artigo).
## Saidas: resultados/calib_controle_result.rds e resultados/resultado_calib_controle.txt
source("copula_ml.R")
sink_file <- file.path(DIR_RESULTADOS, "resultado_calib_controle.txt")
zz <- file(sink_file, open = "wt"); sink(zz, split = TRUE)
R <- readRDS("reais_result_sem.rds")
d <- R$d; X <- R$X; Z <- R$Z; y <- R$y; phaseI <- R$phaseI; n <- nrow(d)
kx <- ncol(X); kz <- ncol(Z)
I <- seq_len(phaseI); II <- (phaseI + 1):n
F_out <- readRDS("farrington_result.rds")$alarms_idx
C_out <- readRDS("cusum_result.rds")$flag
ctrl <- setdiff(II, union(F_out, C_out))
cat("=== CALIBRACAO SOB CONTROLE E DETECCAO NO SURTO (harmonico SEMESTRAL) ===\n")
cat(sprintf("Fase I: %d semanas (%s a %s); Fase II: %d semanas (%s a %s)\n", length(I), d$epiweek[1], d$epiweek[phaseI],
            length(II), d$epiweek[phaseI + 1], d$epiweek[n]))
cat(sprintf("surto: Farrington %d semanas (%s a %s); CUSUM %d; uniao %d; Fase II sob controle: %d semanas\n\n",
            length(F_out), d$epiweek[min(F_out)], d$epiweek[max(F_out)], length(C_out), length(union(F_out, C_out)), length(ctrl)))

fits <- fit_copula_multi_tau(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma",
                             taus = c(0.50, 0.75, 0.90, 0.95, 0.99), verbose = FALSE)
fb <- R$fits[["beta 1 1"]]; fk5 <- R$fits[["kuma 1 1"]]
pred <- function(par, family, tau = 0.5) {
  mu <- plogis(as.numeric(X %*% par[1:kx])); sh <- exp(as.numeric(Z %*% par[(kx + 1):(kx + kz)]))
  ar <- par[kx + kz + 1]; ma <- par[kx + kz + 2]
  Ft <- if (family == "kuma") pkuma(y, mu, sh, tau) else pbetamp(y, mu, sh)
  eps <- qnorm(pmin(pmax(Ft, 1e-12), 1 - 1e-12))
  mom <- .one_step_moments(eps, .dl_recursion(as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = n - 1))))
  Q <- function(u) if (family == "kuma") qkuma(u, mu, sh, tau) else qbetamp(u, mu, sh)
  list(marg = function(t0) Q(rep(t0, n)), cond = function(t0) Q(pnorm(mom$m + sqrt(mom$s2) * qnorm(t0))))
}
pb <- pred(fb$par, "beta")
linha <- function(lim, rot, t0) data.frame(
  limite = rot, tau = t0,
  faseI = mean(y[I] <= lim[I]), faseII_ctrl = mean(y[ctrl] <= lim[ctrl]), faseII_todas = mean(y[II] <= lim[II]),
  exc_faseI = sum(y[I] > lim[I]), exc_ctrl = sum(y[ctrl] > lim[ctrl]), surto_acima = sum(y[F_out] > lim[F_out]),
  lim_medio_I = mean(lim[I]), lim_medio_ctrl = mean(lim[ctrl]))
rows <- list()
for (t0 in c(0.50, 0.75, 0.90, 0.95, 0.99)) {
  pk <- pred(fits[[as.character(t0)]]$par, "kuma", t0)
  rows <- c(rows, list(linha(pk$marg(t0), "Kuma marginal (ajuste em tau)", t0),
                       linha(pb$marg(t0), "Beta marginal (derivado da media)", t0),
                       linha(pk$cond(t0), "Kuma condicional (ajuste em tau)", t0),
                       linha(pb$cond(t0), "Beta condicional (derivado da media)", t0)))
}
## tau = 0.90: segundo ponto estacionario (reinicializacao a partir de tau = 0.95), nota da tab:multitau
nll90 <- function(par) nll_copula(par, y = y[I], X = X[I, ], Z = Z[I, ], p = 1, q = 1, family = "kuma", tau = 0.90)
p95 <- fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma", tau = 0.95)$par
o_de95 <- optim(p95, nll90, method = "BFGS", control = list(maxit = 1500, reltol = 1e-11))
rows <- c(rows, list(linha(pred(o_de95$par, "kuma", 0.90)$marg(0.90), "Kuma marginal, tau=0.90, 2o ponto estacionario", 0.90)))
tab <- do.call(rbind, rows)
cat(sprintf("loglik do 2o ponto estacionario em tau=0.90: %.3f\n", -o_de95$value))
cat(sprintf("alvo de excedencias: Fase I (300) %.1f/%.1f/%.1f; Fase II sob controle (%d) %.1f/%.1f/%.1f para tau = 0.90/0.95/0.99\n\n",
            30, 15, 3, length(ctrl), 0.10 * length(ctrl), 0.05 * length(ctrl), 0.01 * length(ctrl)))
options(width = 200)
print(tab, digits = 3, row.names = FALSE)
## variante sem o termo de tendencia (locacao e forma com intercepto e harmonicos semestrais)
cat("\n--- variante SEM TENDENCIA (limites marginais, tau = 0.90/0.95/0.99) ---\n")
Xs <- cbind(1, d$cos26, d$sin26); Zs <- Xs; ks <- ncol(Xs)
fs <- fit_copula_multi_tau(y[I], Xs[I, ], Zs[I, ], p = 1, q = 1, family = "kuma", taus = c(0.50, 0.75, 0.90, 0.95, 0.99), verbose = FALSE)
fbs <- fit_copula(y[I], Xs[I, ], Zs[I, ], p = 1, q = 1, family = "beta")
mubs <- plogis(as.numeric(Xs %*% fbs$par[1:ks])); phbs <- exp(as.numeric(Zs %*% fbs$par[(ks + 1):(2 * ks)]))
rows_s <- list()
for (t0 in c(0.90, 0.95, 0.99)) {
  rows_s <- c(rows_s, list(linha(plogis(as.numeric(Xs %*% fs[[as.character(t0)]]$par[1:ks])), "Kuma marginal sem tendencia", t0),
                           linha(qbeta(t0, mubs * phbs, (1 - mubs) * phbs), "Beta marginal sem tendencia", t0)))
}
tab_s <- do.call(rbind, rows_s); print(tab_s, digits = 3, row.names = FALSE)
cat(sprintf("AIC (Fase I): Kuma sem tendencia %.2f ; beta sem tendencia %.2f\n", -2 * fs[["0.5"]]$loglik + 2 * length(fs[["0.5"]]$par), fbs$aic))
cat(sprintf("\nquantis empiricos da incidencia na Fase I (0.90/0.95/0.99): %s; na Fase II sob controle: %s\n",
            paste(sprintf("%.4f", quantile(d$incidence[I], c(.90, .95, .99))), collapse = " / "),
            paste(sprintf("%.4f", quantile(d$incidence[ctrl], c(.90, .95, .99))), collapse = " / ")))
cat(sprintf("um caso = 10/3025 = %.5f; maximo da Fase I = %.4f\n", 10 / 3025, max(d$incidence[I])))
sink(); close(zz)
saveRDS(list(tab = tab, tab_sem_tendencia = tab_s, ctrl = ctrl, F_out = F_out, C_out = C_out, loglik_2o = -o_de95$value),
        file.path(DIR_RESULTADOS, "calib_controle_result.rds"))
