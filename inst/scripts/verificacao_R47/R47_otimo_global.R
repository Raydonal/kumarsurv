## R47_otimo_global.R -- o ajuste adotado (Kuma, Z = X, ARMA(1,1), tau = 0.5) e' o maximo global?
## Um submodelo aninhado (mesmo preditor de locacao, forma so' com intercepto) atinge loglik 2543.21
## em analise_harmonicos_inflacao.R, acima dos 2542.27 do modelo completo.
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); d <- R$d; I <- seq_len(R$phaseI); y <- R$y[I]
X <- R$X[I, ]; Z <- R$Z[I, ]; f_adot <- R$fits[["kuma 1 1"]]
obj <- function(p) nll_copula(p, y, X, Z, 1, 1, "kuma", 0.5)
## submodelo: forma so' com intercepto, mesmas covariaveis de locacao
objA <- function(p) nll_copula(c(p[1:4], p[5], 0, 0, 0, p[6:7]), y, X, Z, 1, 1, "kuma", 0.5)
fA <- fit_copula(y, X, matrix(1, length(y), 1), p = 1, q = 1, family = "kuma", tau = 0.5)
cat(sprintf("submodelo (Z = 1): loglik = %.3f  par = %s\n", fA$loglik, paste(sprintf("%.3f", fA$par), collapse = " ")))
cat(sprintf("modelo adotado arquivado: loglik = %.3f\n", f_adot$loglik))
## partida no submodelo, completando os coeficientes de forma com zero
st <- c(fA$par[1:4], fA$par[5], 0, 0, 0, fA$par[6:7])
cat(sprintf("nll do completo na partida (= submodelo): %.3f\n", -obj(st)))
o1 <- optim(st, obj, method = "BFGS", control = list(maxit = 3000, reltol = 1e-12))
o2 <- optim(o1$par, obj, method = "BFGS", control = list(maxit = 3000, reltol = 1e-13))
cat(sprintf("completo reotimizado a partir do submodelo: loglik = %.3f (conv %d)\n", -o2$value, o2$convergence))
print(round(rbind(adotado = f_adot$par, reotimizado = o2$par), 4))
## busca multi-partida simples
set.seed(47); best <- -o2$value
for (k in 1:20) { s0 <- o2$par + rnorm(10, 0, c(1, 3, .5, .5, .3, .5, .1, .1, .02, .05))
  ok <- try(optim(s0, obj, method = "BFGS", control = list(maxit = 3000, reltol = 1e-12)), silent = TRUE)
  if (!inherits(ok, "try-error") && is.finite(ok$value) && -ok$value > best + 1e-6) { best <- -ok$value; cat(sprintf("  partida %d: loglik %.3f\n", k, best)) } }
cat(sprintf("melhor loglik encontrada: %.3f\n", best))
saveRDS(list(fA = fA, o2 = o2), file.path(DIR_RESULTADOS, "R47_otimo_global.rds"))
