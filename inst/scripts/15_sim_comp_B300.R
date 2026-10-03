## 15_sim_comp_B300.R
## Substitui 04_misspec.R (que usava R=50, nao 500 como o .tex afirmava) e
## fornece, pela primeira vez com script proprio, a comparacao de perda
## pinball "estimacao direta em tau" vs "derivada da mediana" vs "beta
## derivada da media" citada em sec:sim-comp mas sem gerador ate esta rodada.
## Mesmo desenho descrito no artigo para os dois experimentos da Secao
## sec:sim-comp: n=250, dependencia ARMA(1,1), B=300 (unificado).
source("copula_ml.R")
n <- 250; set.seed(7)
tt <- ((1:n) - 0.5 * n) / 100
X <- cbind(1, tt, cos(2 * pi * (1:n) / 52), sin(2 * pi * (1:n) / 52))
## Rodada 41: era `Z <- cbind(1)`, matriz 1x1 e nao coluna de n uns (armadilha 8
## da secao 5 do CONTEXTO_RETOMADA). O ajuste sobrevive porque nll_copula e'
## vetorizado e recicla, e por isso a comparacao pinball abaixo nao foi afetada;
## o crps() interno, que indexa sh[i], devolvia NaN — dai o artefato arquivado
## ter sido renomeado _OBSOLETO_CRPS_NaN. O CRPS reportado no artigo vem de
## 17_crps_fix_B300.R; aqui o que interessa e' a perda pinball.
Z <- matrix(1, n, 1)
beta <- c(qlogis(0.45), 0.3, -0.5, -0.2)
ar <- 0.5; ma <- 0.3
gam_k <- log(2.0); gam_b <- log(12)

crps <- function(fit, family, y) {
  par <- fit$par; kx <- ncol(X); kz <- 1
  mu <- plogis(as.numeric(X %*% par[1:kx])); sh <- exp(as.numeric(Z %*% par[(kx + 1):(kx + kz)]))
  arh <- par[kx + kz + 1]; mah <- par[kx + kz + 2]
  if (family == "kuma") { Ft <- pkuma(y, mu, sh); Qf <- function(u, i) qkuma(u, mu[i], sh[i]) }
  else { Ft <- pbetamp(y, mu, sh); Qf <- function(u, i) qbetamp(u, mu[i], sh[i]) }
  eps <- qnorm(pmin(pmax(Ft, 1e-12), 1 - 1e-12))
  rho <- as.numeric(ARMAacf(ar = arh, ma = mah, lag.max = n - 1))
  mom <- .one_step_moments(eps, .dl_recursion(rho))
  m <- mom$m; v <- mom$s2; grid <- seq(0.05, 0.95, 0.05)
  mean(sapply(grid, function(tau) {
    q <- sapply(1:n, function(i) Qf(pnorm(m[i] + sqrt(v[i]) * qnorm(tau)), i))
    mean(ifelse(y >= q, tau * (y - q), (1 - tau) * (q - y)))
  })) * 2
}
pinball <- function(y, q, tau) mean(ifelse(y >= q, tau * (y - q), (1 - tau) * (q - y)))

fitfast <- function(y, fam, tau = 0.5) {
  o0 <- optim(c(qlogis(median(y)), rep(0, ncol(X) - 1 + ncol(Z))),
              function(p) nll_copula(p, y, X, Z, 0, 0, fam, tau = tau),
              method = "BFGS", control = list(maxit = 150))
  o1 <- optim(c(o0$par, 0.01, 0.01),
              function(p) nll_copula(p, y, X, Z, 1, 1, fam, tau = tau),
              method = "BFGS", control = list(maxit = 200, reltol = 1e-9))
  list(par = o1$par, conv = o1$convergence, value = o1$value)
}
## Rodada 47: o ajuste em tau = 0.95/0.99 passa a ter uma segunda partida, com o ARMA do
## ajuste na mediana, e fica com a de maior verossimilhanca. Com partida unica, a replica 20
## terminava, com a biblioteca corrigida, num ponto estacionario 36 unidades de
## log-verossimilhanca abaixo do alcancado pela outra trajetoria (armadilha 29 do
## CONTEXTO_RETOMADA.md). Nas demais replicas as duas partidas coincidem ou a segunda melhora.
fitfast_tau <- function(y, tau, ref) {
  a <- fitfast(y, "kuma", tau = tau)
  if (inherits(ref, "try-error") || ref$conv != 0) return(a)
  o0 <- optim(c(qlogis(median(y)), rep(0, ncol(X) - 1 + ncol(Z))),
              function(p) nll_copula(p, y, X, Z, 0, 0, "kuma", tau = tau),
              method = "BFGS", control = list(maxit = 150))
  b <- optim(c(o0$par, ref$par[ncol(X) + ncol(Z) + 1:2]),
             function(p) nll_copula(p, y, X, Z, 1, 1, "kuma", tau = tau),
             method = "BFGS", control = list(maxit = 200, reltol = 1e-9))
  if (b$convergence == 0 && (a$conv != 0 || b$value < a$value - 1e-6)) list(par = b$par, conv = 0L, value = b$value) else a
}

B <- 300L
res_crps <- array(NA, c(B, 2, 2), dimnames = list(NULL, c("DGP_kuma", "DGP_beta"), c("fit_kuma", "fit_beta")))
pin_direto_95 <- pin_deriv_95 <- pin_beta_95 <- rep(NA, B)
pin_direto_99 <- pin_deriv_99 <- pin_beta_99 <- rep(NA, B)

kx <- ncol(X); kz <- 1
t0 <- Sys.time()
for (i in seq_len(B)) {
  yk <- sim_copula(X, Z, beta, gam_k, ar, ma, family = "kuma")
  yb <- sim_copula(X, Z, beta, gam_b, ar, ma, family = "beta")

  ## --- experimento 1: custo de ma especificacao (CRPS), ambos os DGPs ---
  fk_on_k <- try(suppressWarnings(fitfast(yk, "kuma")), silent = TRUE)
  fb_on_k <- try(suppressWarnings(fitfast(yk, "beta")), silent = TRUE)
  fk_on_b <- try(suppressWarnings(fitfast(yb, "kuma")), silent = TRUE)
  fb_on_b <- try(suppressWarnings(fitfast(yb, "beta")), silent = TRUE)
  if (!inherits(fk_on_k, "try-error") && fk_on_k$conv == 0) res_crps[i, "DGP_kuma", "fit_kuma"] <- crps(fk_on_k, "kuma", yk)
  if (!inherits(fb_on_k, "try-error") && fb_on_k$conv == 0) res_crps[i, "DGP_kuma", "fit_beta"] <- crps(fb_on_k, "beta", yk)
  if (!inherits(fk_on_b, "try-error") && fk_on_b$conv == 0) res_crps[i, "DGP_beta", "fit_kuma"] <- crps(fk_on_b, "kuma", yb)
  if (!inherits(fb_on_b, "try-error") && fb_on_b$conv == 0) res_crps[i, "DGP_beta", "fit_beta"] <- crps(fb_on_b, "beta", yb)

  ## --- experimento 2: estimacao direta em tau vs derivada da mediana/media (DGP Kuma) ---
  ## fk_on_k acima ja e' o ajuste Kuma em tau=0.5 (mediana); fb_on_k ja e' o ajuste beta (media)
  if (!inherits(fk_on_k, "try-error") && fk_on_k$conv == 0) {
    par <- fk_on_k$par; mu_med <- plogis(as.numeric(X %*% par[1:kx])); sh <- exp(as.numeric(Z %*% par[(kx + 1):(kx + kz)]))
    kap <- log(0.5) / log(1 - mu_med^sh)   ## 'a' implicito no ajuste na mediana
    q_deriv_95 <- (1 - (1 - 0.95)^(1 / kap))^(1 / sh)
    q_deriv_99 <- (1 - (1 - 0.99)^(1 / kap))^(1 / sh)
    pin_deriv_95[i] <- pinball(yk, q_deriv_95, 0.95)
    pin_deriv_99[i] <- pinball(yk, q_deriv_99, 0.99)
  }
  if (!inherits(fb_on_k, "try-error") && fb_on_k$conv == 0) {
    par <- fb_on_k$par; mu_b <- plogis(as.numeric(X %*% par[1:kx])); phi_b <- exp(as.numeric(Z %*% par[(kx + 1):(kx + kz)]))
    q_beta_95 <- qbetamp(0.95, mu_b, phi_b); q_beta_99 <- qbetamp(0.99, mu_b, phi_b)
    pin_beta_95[i] <- pinball(yk, q_beta_95, 0.95)
    pin_beta_99[i] <- pinball(yk, q_beta_99, 0.99)
  }
  f95 <- try(suppressWarnings(fitfast_tau(yk, 0.95, fk_on_k)), silent = TRUE)
  f99 <- try(suppressWarnings(fitfast_tau(yk, 0.99, fk_on_k)), silent = TRUE)
  if (!inherits(f95, "try-error") && f95$conv == 0) {
    q_dir_95 <- plogis(as.numeric(X %*% f95$par[1:kx]))
    pin_direto_95[i] <- pinball(yk, q_dir_95, 0.95)
  }
  if (!inherits(f99, "try-error") && f99$conv == 0) {
    q_dir_99 <- plogis(as.numeric(X %*% f99$par[1:kx]))
    pin_direto_99[i] <- pinball(yk, q_dir_99, 0.99)
  }
  if (i %% 25 == 0) cat("i =", i, "/", B, "\n")
}
t1 <- Sys.time()
cat("tempo total:", round(as.numeric(difftime(t1, t0, units = "mins")), 2), "min\n")

M <- apply(res_crps, c(2, 3), mean, na.rm = TRUE)
cat("\n=== CRPS medio (menor e' melhor); linhas=verdade, colunas=modelo ajustado ===\n")
print(round(M, 5))
cat(sprintf("\nDGP Kumaraswamy: beta perde %.1f%% em CRPS\n", 100 * (M["DGP_kuma", "fit_beta"] / M["DGP_kuma", "fit_kuma"] - 1)))
cat(sprintf("DGP Beta:        kuma perde %.1f%% em CRPS\n", 100 * (M["DGP_beta", "fit_kuma"] / M["DGP_beta", "fit_beta"] - 1)))

cat("\n=== Pinball: direto em tau vs derivado da mediana vs beta derivado da media (DGP Kuma) ===\n")
cat(sprintf("tau=0.95: direto=%.5f  derivado(mediana)=%.5f  beta(media)=%.5f  (n_ok=%d/%d/%d)\n",
            mean(pin_direto_95, na.rm = TRUE), mean(pin_deriv_95, na.rm = TRUE), mean(pin_beta_95, na.rm = TRUE),
            sum(!is.na(pin_direto_95)), sum(!is.na(pin_deriv_95)), sum(!is.na(pin_beta_95))))
cat(sprintf("tau=0.99: direto=%.5f  derivado(mediana)=%.5f  beta(media)=%.5f  (n_ok=%d/%d/%d)\n",
            mean(pin_direto_99, na.rm = TRUE), mean(pin_deriv_99, na.rm = TRUE), mean(pin_beta_99, na.rm = TRUE),
            sum(!is.na(pin_direto_99)), sum(!is.na(pin_deriv_99)), sum(!is.na(pin_beta_99))))

## Rodada 47: diferenca pareada direto - derivado, erro-padrao de Monte Carlo e razao t
## (citados na Secao 5.3), e excesso de perda do quantil beta sobre o direto
for (tt in c("95", "99")) {
  d <- get(paste0("pin_direto_", tt)) - get(paste0("pin_deriv_", tt)); ep <- sd(d, na.rm = TRUE) / sqrt(sum(!is.na(d)))
  cat(sprintf("tau=0.%s: diferenca pareada direto - derivado = %.2e, EP de Monte Carlo = %.2e, razao t = %.2f; beta perde %.1f%% em relacao ao direto\n",
              tt, mean(d, na.rm = TRUE), ep, mean(d, na.rm = TRUE) / ep,
              100 * (mean(get(paste0("pin_beta_", tt)), na.rm = TRUE) / mean(get(paste0("pin_direto_", tt)), na.rm = TRUE) - 1)))
}

saveRDS(list(B = B, res_crps = res_crps, M = M,
             pin_direto_95 = pin_direto_95, pin_deriv_95 = pin_deriv_95, pin_beta_95 = pin_beta_95,
             pin_direto_99 = pin_direto_99, pin_deriv_99 = pin_deriv_99, pin_beta_99 = pin_beta_99),
        "sim_comp_B300_result.rds")
