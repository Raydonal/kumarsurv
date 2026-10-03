## 14_mc_validacao_B300.R
## Rodada 47: hessiana numerica com passo relativo d = 1e-3 (o padrao do numDeriv, d = 0.1,
## levava os coeficientes AR para fora da regiao estacionaria e subestimava os EP de dependencia).
## Substitui 05_mc_validacao.R (que usava R=60, sem a coluna "Independente").
## Reproduz a Tabela STK+TAB1: desenho de Varin (2014), ARMA(2,1), n=468,
## harmonico anual, marginal Kumaraswamy, com B=300 réplicas (unificado com os
## demais experimentos do artigo), ajustando tanto o modelo com erros ARMA(2,1)
## quanto o modelo com erros independentes (p=q=0), como na tabela publicada.
source("copula_ml.R")
N <- 52 * 9
time <- ((1:N) - 0.5 * N) / 100
X <- cbind(1, time, cos(2 * pi * (1:N) / 52), sin(2 * pi * (1:N) / 52)); Z <- X
colnames(X) <- colnames(Z) <- c("int", "trend", "cos", "sin")
beta  <- c(-2, 0.5, -0.67, -0.22)
gamma <- c(-2, 0.1, -0.19, -0.06)
ar <- c(1.5, -0.6); ma <- -0.3
tv_arma <- c(beta, gamma, ar, ma)
tv_ind  <- c(beta, gamma)
K_arma <- length(tv_arma); K_ind <- length(tv_ind)

fit_arma <- function(y) {
  o0 <- optim(c(qlogis(median(y)), rep(0, 7)),
              function(p) nll_copula(p, y, X, Z, 0, 0, "kuma"),
              method = "BFGS", control = list(maxit = 200))
  o1 <- optim(c(o0$par, rep(0.01, 3)),
              function(p) nll_copula(p, y, X, Z, 2, 1, "kuma"),
              method = "BFGS", control = list(maxit = 300, reltol = 1e-9))
  list(par = o1$par, conv = o1$convergence,
       se = tryCatch(sqrt(diag(solve(numDeriv::hessian(
         function(p) nll_copula(p, y, X, Z, 2, 1, "kuma"), o1$par, method.args = list(d = 1e-3))))),
         error = function(e) rep(NA, length(o1$par))))
}
fit_ind <- function(y) {
  o0 <- optim(c(qlogis(median(y)), rep(0, 7)),
              function(p) nll_copula(p, y, X, Z, 0, 0, "kuma"),
              method = "BFGS", control = list(maxit = 300, reltol = 1e-10))
  list(par = o0$par, conv = o0$convergence,
       se = tryCatch(sqrt(diag(solve(numDeriv::hessian(
         function(p) nll_copula(p, y, X, Z, 0, 0, "kuma"), o0$par, method.args = list(d = 1e-3))))),
         error = function(e) rep(NA, length(o0$par))))
}

B <- 300L
est_arma <- matrix(NA, B, K_arma); se_arma <- matrix(NA, B, K_arma)
est_ind  <- matrix(NA, B, K_ind);  se_ind  <- matrix(NA, B, K_ind)
set.seed(2024)
conv_arma <- 0L; conv_ind <- 0L
t0 <- Sys.time()
for (i in seq_len(B)) {
  y <- sim_copula(X, Z, beta, gamma, ar, ma, family = "kuma")
  fa <- try(suppressWarnings(fit_arma(y)), silent = TRUE)
  if (!inherits(fa, "try-error") && fa$conv == 0) {
    est_arma[i, ] <- fa$par; se_arma[i, ] <- fa$se; conv_arma <- conv_arma + 1L
  }
  fi <- try(suppressWarnings(fit_ind(y)), silent = TRUE)
  if (!inherits(fi, "try-error") && fi$conv == 0) {
    est_ind[i, ] <- fi$par; se_ind[i, ] <- fi$se; conv_ind <- conv_ind + 1L
  }
  if (i %% 25 == 0) cat("i =", i, "/", B,
                        " conv_arma =", conv_arma, " conv_ind =", conv_ind, "\n")
}
t1 <- Sys.time()
cat("tempo total:", round(as.numeric(difftime(t1, t0, units = "mins")), 2), "min\n")

nm_arma <- c(paste0("mu.", colnames(X)), paste0("sh.", colnames(Z)), "ar1", "ar2", "ma1")
nm_ind  <- c(paste0("mu.", colnames(X)), paste0("sh.", colnames(Z)))

## Rodada 47: desvio-padrao robusto (amplitude interquartil / 1.349), coluna DPr do
## Quadro 8, e as matrizes de EP por replica passam a ser gravadas.
dprob <- function(E) apply(E, 2, function(x) IQR(x, na.rm = TRUE) / 1.349)
tab <- data.frame(
  param = nm_arma,
  verdadeiro = round(tv_arma, 3),
  media_arma = round(colMeans(est_arma, na.rm = TRUE), 3),
  dp_arma = round(apply(est_arma, 2, sd, na.rm = TRUE), 3),
  dprob_arma = round(dprob(est_arma), 3),
  ep_arma = round(colMeans(se_arma, na.rm = TRUE), 3),
  media_ind = c(round(colMeans(est_ind, na.rm = TRUE), 3), NA, NA, NA),
  dp_ind = c(round(apply(est_ind, 2, sd, na.rm = TRUE), 3), NA, NA, NA),
  dprob_ind = c(round(dprob(est_ind), 3), NA, NA, NA),
  ep_ind = c(round(colMeans(se_ind, na.rm = TRUE), 3), NA, NA, NA)
)
cat("\n=== Tabela STK+TAB1 (B=300) ===\n")
print(tab, row.names = FALSE)
cat("\nConvergidas: ARMA(2,1) =", conv_arma, "/", B, "; Independente =", conv_ind, "/", B, "\n")
cat("EP calculados: ARMA(2,1)", paste(colSums(is.finite(se_arma)), collapse = " "), "\n")
cat("razao DPr/EP, ARMA(2,1):", sprintf("%.2f", tab$dprob_arma / tab$ep_arma), "\n")
cat("razao DPr/EP, independente:", sprintf("%.2f", tab$dprob_ind[1:K_ind] / tab$ep_ind[1:K_ind]), "\n")
ok <- stats::complete.cases(est_arma); Eok <- est_arma[ok, ]
atip <- rowSums(abs(sweep(Eok, 2, apply(Eok, 2, median))) >
                5 * matrix(dprob(Eok), nrow(Eok), ncol(Eok), byrow = TRUE)) > 0
cat(sprintf("replicas com alguma estimativa a mais de 5 DPr da mediana: %d de %d\n", sum(atip), nrow(Eok)))
cat("DP sem essas replicas:", sprintf("%.3f", apply(Eok[!atip, , drop = FALSE], 2, sd)), "\n")
saveRDS(list(tab = tab, B = B, conv_arma = conv_arma, conv_ind = conv_ind,
             est_arma = est_arma, est_ind = est_ind, se_arma = se_arma, se_ind = se_ind),
        "mc_validacao_B300_result.rds")
