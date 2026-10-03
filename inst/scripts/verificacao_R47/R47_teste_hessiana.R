## R47_teste_hessiana.R -- efeito do passo da hessiana numerica (numDeriv, d = 0.1 por
## padrao) sobre os erros-padrao, na aplicacao e numa replica do desenho de Varin.
source("copula_ml.R")
se_of <- function(H) { s <- tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA, nrow(H))); s[!is.finite(s)] <- NA; s }
R <- readRDS("reais_result_sem.rds"); f <- R$fits[["kuma 1 1"]]; I <- seq_len(R$phaseI)
obj <- function(p) nll_copula(p, R$y[I], R$X[I, ], R$Z[I, ], 1, 1, "kuma", 0.5)
cat("aplicacao, parametros:", sprintf("%.4f", f$par), "\n")
cat("penalidade atingida no passo padrao? nll(ar*1.1) =", obj(replace(f$par, 9, f$par[9] * 1.1)), "\n")
H0 <- numDeriv::hessian(obj, f$par)                                   # como no pacote (d = 0.1)
H1 <- numDeriv::hessian(obj, f$par, method.args = list(d = 1e-3))    # passo relativo 0.1%
H2 <- optimHess(f$par, obj, control = list(ndeps = rep(1e-4, length(f$par))))
tab <- rbind(pacote_d0.1 = se_of(H0), numDeriv_d1e.3 = se_of(H1), optimHess_1e.4 = se_of(H2))
colnames(tab) <- names(f$par); print(round(tab, 4))
cat("EP arquivado (fit$se):", sprintf("%.4f", f$se), "\n")
## replica do desenho de Varin
N <- 52 * 9; time <- ((1:N) - 0.5 * N) / 100
X <- cbind(1, time, cos(2 * pi * (1:N) / 52), sin(2 * pi * (1:N) / 52)); Z <- X
tv <- c(-2, 0.5, -0.67, -0.22, -2, 0.1, -0.19, -0.06, 1.5, -0.6, -0.3)
set.seed(2024)
y <- sim_copula(X, Z, tv[1:4], tv[5:8], ar = tv[9:10], ma = tv[11], family = "kuma")
o0 <- optim(c(qlogis(median(y)), rep(0, 7)), function(p) nll_copula(p, y, X, Z, 0, 0, "kuma"), method = "BFGS", control = list(maxit = 200))
ob <- function(p) nll_copula(p, y, X, Z, 2, 1, "kuma")
o1 <- optim(c(o0$par, rep(0.01, 3)), ob, method = "BFGS", control = list(maxit = 300, reltol = 1e-9))
cat("\nVarin (1 replica), estimativas:", sprintf("%.3f", o1$par), "\n")
tb <- rbind(pacote_d0.1 = se_of(numDeriv::hessian(ob, o1$par)), numDeriv_d1e.3 = se_of(numDeriv::hessian(ob, o1$par, method.args = list(d = 1e-3))))
print(round(tb, 4))
cat("DP robusto do Monte Carlo arquivado (ar1, ar2, ma1): 0.0964 0.0862 0.1183; EP medio publicado: 0.017 0.007 0.045\n")
