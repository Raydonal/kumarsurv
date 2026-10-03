## R47_vero_exata.R -- compara nll_copula com a log-verossimilhanca exata (Durbin-Levinson
## completo, sem truncamento), perto da fronteira de estacionariedade. Escrito durante a
## auditoria, com a biblioteca anterior (variancia por 2000 pesos psi, log-determinante com
## L = 80), onde as duas divergiam junto a fronteira; com a biblioteca corrigida da Rodada
## 47 elas coincidem. A funcao kal() abaixo reconstroi de proposito o calculo anterior.
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); I <- seq_len(R$phaseI); y <- R$y[I]; X <- R$X[I, ]; Z <- R$Z[I, ]
f <- R$fits[["kuma 1 1"]]; P <- readRDS(file.path(DIR_RESULTADOS, "R47_perfil_ar1_global.rds"))
ll_exata <- function(par) {
  mu <- plogis(as.numeric(X %*% par[1:4])); sh <- exp(as.numeric(Z %*% par[5:8]))
  eps <- qnorm(pmin(pmax(pkuma(y, mu, sh, 0.5), 1e-12), 1 - 1e-12)); n <- length(y)
  rho <- as.numeric(ARMAacf(ar = par[9], ma = par[10], lag.max = n - 1))
  mom <- .one_step_moments(eps, .dl_recursion(rho))
  sum(dkuma(y, mu, sh, 0.5, log = TRUE)) - 0.5 * (sum(log(mom$s2)) + sum((eps - mom$m)^2 / mom$s2) - sum(eps^2))
}
ll_pacote <- function(par) -nll_copula(par, y, X, Z, 1, 1, "kuma", 0.5)
kal <- function(par, ss) { mu <- plogis(as.numeric(X %*% par[1:4])); sh <- exp(as.numeric(Z %*% par[5:8]))
  eps <- qnorm(pmin(pmax(pkuma(y, mu, sh, 0.5), 1e-12), 1 - 1e-12)); n <- length(y)
  psi <- ARMAtoMA(ar = par[9], ma = par[10], lag.max = 2000); g0 <- 1 + sum(psi^2)
  kr <- KalmanRun(eps * sqrt(g0), makeARIMA(phi = par[9], theta = par[10], Delta = numeric(0), SSinit = ss), update = FALSE)
  n * kr$values[2] }
pts <- list("adotado (0.948)" = f$par, "melhor (0.9995)" = P$best$par)
for (nm in names(pts)) { p <- pts[[nm]]
  cat(sprintf("%-18s pacote: %.3f | exata: %.3f | forma quad. Kalman Gardner %.3f Rossignol %.3f\n", nm, ll_pacote(p), ll_exata(p),
              kal(p, "Gardner1980"), kal(p, "Rossignol2011"))) }
## soma infinita de psi^2 truncada em 2000 perto da raiz unitaria?
p <- P$best$par; psi <- ARMAtoMA(ar = p[9], ma = p[10], lag.max = 2000)
cat(sprintf("\npsi_2000 = %.3e ; variancia exata do ARMA(1,1) = %.4f ; 1 + sum(psi^2) truncada = %.4f\n",
            psi[2000], (1 + 2 * p[9] * p[10] + p[10]^2) / (1 - p[9]^2), 1 + sum(psi^2)))
cat("\nerro do pacote ao variar so' o AR (demais parametros no ajuste adotado):\n")
for (a in c(0.90, 0.95, 0.97, 0.98, 0.985, 0.99, 0.993, 0.996, 0.999)) { p <- f$par; p[9] <- a
  cat(sprintf("  AR=%.3f  pacote %.3f  exata %.3f  erro %+.4f\n", a, ll_pacote(p), ll_exata(p), ll_pacote(p) - ll_exata(p))) }
## submodelo da comparacao de periodicidade (forma so' com intercepto): valor exato no seu otimo
G <- readRDS(file.path(DIR_RESULTADOS, "R47_otimo_global.rds")); pA <- G$fA$par
pfull <- c(pA[1:4], pA[5], 0, 0, 0, pA[6:7])
cat(sprintf("\nsubmodelo semestral (Z = 1), AR = %.4f: pacote %.3f  exata %.3f\n", pA[6], ll_pacote(pfull), ll_exata(pfull)))
