## R47_ep_varin.R -- Rodada 47. Recalcula os erros-padrao assintoticos do experimento
## do desenho de Varin (14_mc_validacao_B300.R) nas estimativas ARQUIVADAS, com passo
## de hessiana numerica adequado (numDeriv, d = 1e-3). O script 14 usava o passo padrao
## do numDeriv (d = 0.1, isto e, 10% do valor do parametro), que leva os coeficientes
## AR para fora da regiao estacionaria, onde nll_copula devolve a penalidade 1e10.
## Os dados de cada replica sao regenerados na mesma sequencia (set.seed(2024) e
## sim_copula no laco, unica fonte de numeros aleatorios); VERIFICA=1 reajusta as
## primeiras replicas e confere com as estimativas arquivadas.
source("copula_ml.R")
m <- readRDS("mc_validacao_B300_result.rds")
N <- 52 * 9; time <- ((1:N) - 0.5 * N) / 100
X <- cbind(1, time, cos(2 * pi * (1:N) / 52), sin(2 * pi * (1:N) / 52)); Z <- X
beta <- c(-2, 0.5, -0.67, -0.22); gamma <- c(-2, 0.1, -0.19, -0.06); ar <- c(1.5, -0.6); ma <- -0.3
i0 <- as.integer(Sys.getenv("I0", "1")); i1 <- as.integer(Sys.getenv("I1", "300"))
out <- file.path(DIR_RESULTADOS, sprintf("R47_ep_varin_%03d_%03d.rds", i0, i1))
set.seed(2024)
se_a <- matrix(NA, 300, 11); se_i <- matrix(NA, 300, 8)
for (i in seq_len(i1)) {
  y <- sim_copula(X, Z, beta, gamma, ar, ma, family = "kuma")
  if (i < i0) next
  if (Sys.getenv("VERIFICA") == "1" && i <= i0 + 1) {
    o0 <- optim(c(qlogis(median(y)), rep(0, 7)), function(p) nll_copula(p, y, X, Z, 0, 0, "kuma"), method = "BFGS", control = list(maxit = 200))
    o1 <- optim(c(o0$par, rep(0.01, 3)), function(p) nll_copula(p, y, X, Z, 2, 1, "kuma"), method = "BFGS", control = list(maxit = 300, reltol = 1e-9))
    cat(sprintf("replica %d: max |reajuste - arquivado| = %.2e\n", i, max(abs(o1$par - m$est_arma[i, ]))))
  }
  if (all(is.finite(m$est_arma[i, ])))
    se_a[i, ] <- tryCatch(sqrt(diag(solve(numDeriv::hessian(function(p) nll_copula(p, y, X, Z, 2, 1, "kuma"),
                           m$est_arma[i, ], method.args = list(d = 1e-3))))), error = function(e) rep(NA, 11))
  if (all(is.finite(m$est_ind[i, ])))
    se_i[i, ] <- tryCatch(sqrt(diag(solve(numDeriv::hessian(function(p) nll_copula(p, y, X, Z, 0, 0, "kuma"),
                           m$est_ind[i, ], method.args = list(d = 1e-3))))), error = function(e) rep(NA, 8))
}
saveRDS(list(i0 = i0, i1 = i1, se_arma = se_a[i0:i1, , drop = FALSE], se_ind = se_i[i0:i1, , drop = FALSE]), out)
cat("gravado", out, "\n")
