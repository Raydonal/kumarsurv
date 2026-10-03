## R47_global_exata.R -- maximo global e perfil do AR(1) com a log-verossimilhanca EXATA
## (Durbin-Levinson completo), restrita a |AR| <= 0.997 e |MA| < 1.
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); I <- seq_len(R$phaseI); y <- R$y[I]; X <- R$X[I, ]; Z <- R$Z[I, ]
f <- R$fits[["kuma 1 1"]]; n <- length(y)
nll_ex <- function(par, fam = "kuma") {
  if (abs(par[9]) > 0.997 || abs(par[10]) >= 1) return(1e10)
  mu <- plogis(as.numeric(X %*% par[1:4])); sh <- exp(as.numeric(Z %*% par[5:8]))
  if (any(!is.finite(mu)) || any(mu <= 0 | mu >= 1) || any(!is.finite(sh))) return(1e10)
  if (fam == "kuma") { ld <- dkuma(y, mu, sh, 0.5, log = TRUE); u <- pkuma(y, mu, sh, 0.5) } else { ld <- dbetamp(y, mu, sh, log = TRUE); u <- pbetamp(y, mu, sh) }
  eps <- qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
  mom <- try(.one_step_moments(eps, .dl_recursion(as.numeric(ARMAacf(ar = par[9], ma = par[10], lag.max = n - 1)))), silent = TRUE)
  if (inherits(mom, "try-error") || any(!is.finite(mom$s2)) || any(mom$s2 <= 0)) return(1e10)
  v <- -(sum(ld) - 0.5 * (sum(log(mom$s2)) + sum((eps - mom$m)^2 / mom$s2) - sum(eps^2)))
  if (!is.finite(v)) 1e10 else v
}
for (fam in c("kuma", "beta")) {
  f0 <- R$fits[[paste(fam, "1 1")]]
  set.seed(4747); st <- list(f0$par)
  for (k in 1:30) st[[length(st) + 1]] <- f0$par + rnorm(10, 0, c(1.5, 4, .5, .5, .4, .8, .1, .1, .02, .05))
  G <- readRDS(file.path(DIR_RESULTADOS, "R47_otimo_global.rds"))
  if (fam == "kuma") st[[length(st) + 1]] <- pmin(pmax(G$o2$par, -50), 50)
  o <- lapply(st, function(s) try(optim(s, nll_ex, fam = fam, method = "BFGS", control = list(maxit = 3000, reltol = 1e-11)), silent = TRUE))
  o <- o[!sapply(o, inherits, "try-error")]; ll <- sapply(o, function(z) -z$value)
  b <- o[[which.max(ll)]]
  cat(sprintf("\n[%s] ajuste arquivado: loglik exata %.3f (AR %.4f)\n", fam, -nll_ex(f0$par, fam), f0$par[9]))
  cat(sprintf("[%s] melhor com verossimilhanca exata: loglik %.3f\n", fam, -b$value)); print(round(b$par, 4))
  cat(sprintf("[%s] maximos distintos (loglik, AR): %s\n", fam, paste(unique(sprintf("%.2f/%.3f", round(ll, 2), sapply(o, function(z) z$par[9])))[order(-ll)][1:8], collapse = "  ")))
  assign(paste0("best_", fam), b)
}
## perfil exato do AR (Kuma), partidas nas duas bacias
pf <- function(phi0) { bl <- -Inf
  for (s in list(f$par[-9], best_kuma$par[-9])) { nf <- function(q) { p <- numeric(10); p[-9] <- q; p[9] <- phi0; nll_ex(p) }
    o <- try(optim(s, nf, method = "BFGS", control = list(maxit = 3000, reltol = 1e-11)), silent = TRUE)
    if (!inherits(o, "try-error") && -o$value > bl) bl <- -o$value }
  bl }
grid <- c(seq(0.85, 0.98, by = 0.01), 0.985, 0.99, 0.993, 0.995, 0.997)
lp <- sapply(grid, pf); lmax <- max(lp, -best_kuma$value)
cat("\nperfil exato do AR (Kuma):\n"); print(round(cbind(AR = grid, loglik = lp, dev = lmax - lp), 3))
cat(sprintf("regiao de 95%% (dev <= 1.921): %s\n", paste(grid[lmax - lp <= 1.921], collapse = " ")))
saveRDS(list(best_kuma = best_kuma, best_beta = best_beta, grid = grid, lp = lp), file.path(DIR_RESULTADOS, "R47_global_exata.rds"))
