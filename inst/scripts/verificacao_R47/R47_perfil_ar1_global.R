## R47_perfil_ar1_global.R -- perfil de verossimilhanca do AR(1) da copula com partidas
## nas duas bacias (a do ajuste adotado e a de maior verossimilhanca, perto da fronteira),
## grade estendida ate 0.999 e restricao |MA| < 1.
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); I <- seq_len(R$phaseI); y <- R$y[I]; X <- R$X[I, ]; Z <- R$Z[I, ]
f <- R$fits[["kuma 1 1"]]; G <- readRDS(file.path(DIR_RESULTADOS, "R47_otimo_global.rds"))
obj <- function(p) nll_copula(p, y, X, Z, 1, 1, "kuma", 0.5)
## melhor ajuste global com varias partidas
set.seed(470); starts <- list(f$par, G$o2$par)
for (k in 1:40) starts[[length(starts) + 1]] <- (if (k %% 2) f$par else G$o2$par) + rnorm(10, 0, c(1.5, 4, .5, .5, .4, .8, .1, .1, .03, .06))
fits <- lapply(starts, function(s) try(optim(s, obj, method = "BFGS", control = list(maxit = 4000, reltol = 1e-12)), silent = TRUE))
fits <- fits[!sapply(fits, inherits, "try-error")]; ll <- sapply(fits, function(o) -o$value)
ok <- sapply(fits, function(o) abs(o$par[10]) < 1 && abs(o$par[9]) < 1)
best <- fits[[which(ok)[which.max(ll[ok])]]]
cat(sprintf("melhor ajuste (|AR|,|MA|<1): loglik %.3f\n", -best$value)); print(round(best$par, 4))
lls <- sort(unique(round(ll[ok], 2)), decreasing = TRUE); cat("maximos locais distintos (loglik):", head(lls, 8), "\n")
## perfil
perf <- function(phi0) {
  nf <- function(pf) { p <- numeric(10); p[-9] <- pf; p[9] <- phi0; obj(p) }
  bl <- -Inf; bma <- NA
  for (s in list(f$par[-9], best$par[-9], G$o2$par[-9])) {
    o <- try(optim(s, nf, method = "BFGS", control = list(maxit = 3000, reltol = 1e-11)), silent = TRUE)
    if (!inherits(o, "try-error") && abs(o$par[9]) < 1 && -o$value > bl) { bl <- -o$value; bma <- o$par[9] }
  }
  c(phi = phi0, ll = bl, ma = bma)
}
grid <- c(seq(0.88, 0.98, by = 0.01), 0.985, 0.99, 0.993, 0.996, 0.999)
P <- t(sapply(grid, perf)); lmax <- max(c(P[, "ll"], -best$value), na.rm = TRUE)
P <- cbind(P, dev = lmax - P[, "ll"]); print(round(P, 3))
cat(sprintf("\nregiao de 95%% (dev <= 1.921): %s\n", paste(P[P[, "dev"] <= 1.921, "phi"], collapse = " ")))
saveRDS(list(best = best, perfil = P, maximos = lls), file.path(DIR_RESULTADOS, "R47_perfil_ar1_global.rds"))
