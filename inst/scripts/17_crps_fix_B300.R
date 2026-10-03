## 17_crps_fix_B300.R
## Corrige um bug real encontrado em 04_misspec.R (herdado, presente tambem em
## 15_sim_comp_B300.R por copia fiel): `Z<-cbind(1)` cria matriz 1x1 em R, nao
## uma coluna de 1's de comprimento n. Isso nao afeta o ajuste (nll_copula usa
## sh em aritmetica vetorizada, que recicla length-1 corretamente), mas quebra
## crps()/Qf, que indexa sh[i] explicitamente — sh[i] para i>1 retorna NA.
## Correcao: Z <- matrix(1, n, 1). Recalcula SOMENTE a matriz de CRPS (a
## comparacao pinball de 15_sim_comp_B300.R nao foi afetada, pois usa apenas
## aritmetica vetorizada, sem indexacao sh[i]).
source("copula_ml.R")
n <- 250; set.seed(7)
tt <- ((1:n) - 0.5 * n) / 100
X <- cbind(1, tt, cos(2 * pi * (1:n) / 52), sin(2 * pi * (1:n) / 52))
Z <- matrix(1, n, 1)  ## CORRIGIDO (era cbind(1))
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
fitfast <- function(y, fam) {
  o0 <- optim(c(qlogis(median(y)), rep(0, ncol(X) - 1 + ncol(Z))),
              function(p) nll_copula(p, y, X, Z, 0, 0, fam),
              method = "BFGS", control = list(maxit = 150))
  o1 <- optim(c(o0$par, 0.01, 0.01),
              function(p) nll_copula(p, y, X, Z, 1, 1, fam),
              method = "BFGS", control = list(maxit = 200, reltol = 1e-9))
  list(par = o1$par, conv = o1$convergence)
}

B <- 300L
res_crps <- array(NA, c(B, 2, 2), dimnames = list(NULL, c("DGP_kuma", "DGP_beta"), c("fit_kuma", "fit_beta")))
t0 <- Sys.time()
for (i in seq_len(B)) {
  yk <- sim_copula(X, Z, beta, gam_k, ar, ma, family = "kuma")
  yb <- sim_copula(X, Z, beta, gam_b, ar, ma, family = "beta")
  fk_on_k <- try(suppressWarnings(fitfast(yk, "kuma")), silent = TRUE)
  fb_on_k <- try(suppressWarnings(fitfast(yk, "beta")), silent = TRUE)
  fk_on_b <- try(suppressWarnings(fitfast(yb, "kuma")), silent = TRUE)
  fb_on_b <- try(suppressWarnings(fitfast(yb, "beta")), silent = TRUE)
  if (!inherits(fk_on_k, "try-error") && fk_on_k$conv == 0) res_crps[i, "DGP_kuma", "fit_kuma"] <- crps(fk_on_k, "kuma", yk)
  if (!inherits(fb_on_k, "try-error") && fb_on_k$conv == 0) res_crps[i, "DGP_kuma", "fit_beta"] <- crps(fb_on_k, "beta", yk)
  if (!inherits(fk_on_b, "try-error") && fk_on_b$conv == 0) res_crps[i, "DGP_beta", "fit_kuma"] <- crps(fk_on_b, "kuma", yb)
  if (!inherits(fb_on_b, "try-error") && fb_on_b$conv == 0) res_crps[i, "DGP_beta", "fit_beta"] <- crps(fb_on_b, "beta", yb)
  if (i %% 25 == 0) cat("i =", i, "/", B, "\n")
}
t1 <- Sys.time()
cat("tempo total:", round(as.numeric(difftime(t1, t0, units = "mins")), 2), "min\n")
nNA <- apply(res_crps, c(2, 3), function(v) sum(is.na(v)))
cat("\nNA por celula (de", B, "):\n"); print(nNA)
M <- apply(res_crps, c(2, 3), mean, na.rm = TRUE)
cat("\n=== CRPS medio (menor e' melhor); linhas=verdade, colunas=modelo ajustado ===\n")
print(round(M, 5))
cat(sprintf("\nDGP Kumaraswamy: beta perde %.2f%% em CRPS\n", 100 * (M["DGP_kuma", "fit_beta"] / M["DGP_kuma", "fit_kuma"] - 1)))
cat(sprintf("DGP Beta:        kuma perde %.2f%% em CRPS\n", 100 * (M["DGP_beta", "fit_kuma"] / M["DGP_beta", "fit_beta"] - 1)))
saveRDS(list(B = B, res_crps = res_crps, M = M), "crps_fix_B300_result.rds")
