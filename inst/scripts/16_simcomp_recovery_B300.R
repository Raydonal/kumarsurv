## 16_simcomp_recovery_B300.R
## Substitui 12_simulacao_comparativa.R (R=120) por B=300 (unificado com os
## demais experimentos). Recuperacao de parametros sob cada marginal
## (Kumaraswamy e beta), mesmo desenho, cada modelo ajustado a dados gerados
## pela sua propria marginal. Gera a Tabela tab:simrec.
source("copula_ml.R")
n <- 250; tt <- ((1:n) - n / 2) / 100
X <- cbind(1, tt, cos(2 * pi * (1:n) / 52), sin(2 * pi * (1:n) / 52)); Z <- matrix(1, n, 1)
colnames(X) <- c("int", "trend", "cos", "sin"); colnames(Z) <- "int"
beta <- c(qlogis(0.30), 0.25, -0.45, -0.18)
gk <- log(2.0); gb <- log(15)
ar <- 0.6; ma <- 0.3
fitq <- function(y, fam) {
  o0 <- optim(c(qlogis(median(y)), rep(0, 4)),
              function(p) nll_copula(p, y, X, Z, 0, 0, fam),
              method = "BFGS", control = list(maxit = 300))
  optim(c(o0$par, 0.01, 0.01),
        function(p) nll_copula(p, y, X, Z, 1, 1, fam),
        method = "BFGS", control = list(maxit = 400, reltol = 1e-9))
}
B <- 300L
run <- function(fam, gtrue) {
  tv <- c(beta, gtrue, ar, ma); K <- length(tv); est <- matrix(NA, B, K)
  set.seed(99)
  t0 <- Sys.time()
  for (i in seq_len(B)) {
    y <- sim_copula(X, Z, beta, gtrue, ar, ma, family = fam)
    f <- try(suppressWarnings(fitq(y, fam)), silent = TRUE)
    if (!inherits(f, "try-error") && f$convergence == 0) est[i, ] <- f$par
    if (i %% 50 == 0) cat("  [", fam, "] i =", i, "/", B, "\n")
  }
  t1 <- Sys.time()
  cat("  tempo", fam, ":", round(as.numeric(difftime(t1, t0, units = "mins")), 2), "min\n")
  mn <- colMeans(est, na.rm = TRUE); sdv <- apply(est, 2, sd, na.rm = TRUE)
  data.frame(par = c("mu:int", "mu:trend", "mu:cos", "mu:sin", "forma", "ar1", "ma1"),
             verdadeiro = round(tv, 3), media = round(mn, 3),
             vies = round(mn - tv, 3), EQM = round((mn - tv)^2 + sdv^2, 4),
             conv = sum(!is.na(est[, 1])))
}
cat("=== Kumaraswamy (marginal correta), B=300 ===\n"); rk <- run("kuma", gk); print(rk, row.names = FALSE)
cat("\n=== Beta (marginal correta), B=300 ===\n");        rb <- run("beta", gb); print(rb, row.names = FALSE)
saveRDS(list(kuma = rk, beta = rb, B = B), "simcomp_recovery_B300_result.rds")
