## 13_bootstrap_gof.R
## Bootstrap parametrico (Algoritmo 1 do artigo) para as estatisticas de
## Cramer-von Mises (W^2) e Anderson-Darling (A^2) sobre a PIT dos residuos
## quantilicos preditivos, para as marginais Kumaraswamy e beta ajustadas na
## Fase I (n=300, harmonico semestral, ARMA(1,1)). Nunca executado em nenhuma
## rodada anterior: os valores W^2/A^2 na Tabela tab:comp vinham do documento
## original dos autores, sem p-valor de bootstrap associado.
##
## Uso: B replicas simuladas do modelo AJUSTADO (theta_hat fixo), cada uma
## reestimada por maxima verossimilhanca (theta_hat^(b)), reproduzindo o
## mesmo efeito de nao-pivotalidade presente na estatistica observada.
## p-valor = (1 + #{T_b >= T_obs}) / (B + 1).

source("copula_ml.R")
set.seed(2026)

d <- read.csv("dda_platina_sem.csv", stringsAsFactors = FALSE)
y_raw <- d$incidence
y <- ifelse(y_raw == 0, 1e-7, y_raw)
X <- cbind(1, d$trend, d$cos26, d$sin26); Z <- X
colnames(X) <- colnames(Z) <- c("int", "trend", "cos26", "sin26")
idx <- 1:300
yI <- y[idx]; XI <- X[idx, ]; ZI <- Z[idx, ]
n <- length(yI)

B <- 300L  ## unificado com os demais experimentos de simulacao do artigo

## --- estatisticas EDF sobre a PIT u_(1) <= ... <= u_(n) (formulas do artigo) ---
cvm_stat <- function(u) {
  u <- sort(u); i <- seq_len(length(u))
  1 / (12 * length(u)) + sum((u - (2 * i - 1) / (2 * length(u)))^2)
}
ad_stat <- function(u) {
  n <- length(u); us <- sort(u); i <- seq_len(n)
  us <- pmin(pmax(us, 1e-12), 1 - 1e-12)
  -n - (1 / n) * sum((2 * i - 1) * (log(us) + log(1 - rev(us))))
}

run_family <- function(family) {
  cat("\n=== Familia:", family, "===\n")
  t0 <- Sys.time()
  f0 <- fit_copula(yI, XI, ZI, p = 1, q = 1, family = family)
  qr0 <- qresiduals(f0)
  u0 <- pnorm(qr0$r)
  W2_obs <- cvm_stat(u0); A2_obs <- ad_stat(u0)
  cat("Ajuste observado: loglik =", round(f0$loglik, 3),
      " conv =", f0$convergence,
      " W2_obs =", round(W2_obs, 4), " A2_obs =", round(A2_obs, 4), "\n")

  kx <- ncol(XI); kz <- ncol(ZI)
  beta_hat  <- f0$par[1:kx]
  gamma_hat <- f0$par[(kx + 1):(kx + kz)]
  ar_hat <- f0$par[kx + kz + 1]
  ma_hat <- f0$par[kx + kz + 2]

  W2_b <- numeric(B); A2_b <- numeric(B); conv_b <- logical(B)
  for (b in seq_len(B)) {
    yb <- try(sim_copula(XI, ZI, beta_hat, gamma_hat, ar = ar_hat, ma = ma_hat,
                          family = family), silent = TRUE)
    if (inherits(yb, "try-error") || any(!is.finite(yb)) || any(yb <= 0 | yb >= 1)) {
      conv_b[b] <- FALSE; W2_b[b] <- NA; A2_b[b] <- NA; next
    }
    fb <- try(fit_copula(yb, XI, ZI, p = 1, q = 1, family = family), silent = TRUE)
    if (inherits(fb, "try-error") || fb$convergence != 0) {
      conv_b[b] <- FALSE; W2_b[b] <- NA; A2_b[b] <- NA; next
    }
    qrb <- qresiduals(fb); ub <- pnorm(qrb$r)
    W2_b[b] <- cvm_stat(ub); A2_b[b] <- ad_stat(ub); conv_b[b] <- TRUE
    if (b %% 25 == 0) cat("  b =", b, "/", B,
                          " (convergidas ate agora:", sum(conv_b[1:b]), ")\n")
  }
  t1 <- Sys.time()
  ok <- conv_b
  pW2 <- (1 + sum(W2_b[ok] >= W2_obs, na.rm = TRUE)) / (sum(ok) + 1)
  pA2 <- (1 + sum(A2_b[ok] >= A2_obs, na.rm = TRUE)) / (sum(ok) + 1)
  cat("Convergidas:", sum(ok), "/", B, "\n")
  cat("p-valor bootstrap W2 =", round(pW2, 4), "\n")
  cat("p-valor bootstrap A2 =", round(pA2, 4), "\n")
  cat("tempo total:", round(as.numeric(difftime(t1, t0, units = "mins")), 2), "min\n")
  list(family = family, fit0 = f0, W2_obs = W2_obs, A2_obs = A2_obs,
       W2_b = W2_b, A2_b = A2_b, conv_b = conv_b, B = B,
       pW2 = pW2, pA2 = pA2)
}

res_kuma <- run_family("kuma")
res_beta <- run_family("beta")

saveRDS(list(kuma = res_kuma, beta = res_beta, B = B),
        "bootstrap_gof_result.rds")

cat("\n=== RESUMO FINAL ===\n")
cat(sprintf("Kumaraswamy: W2_obs=%.4f (p=%.4f)  A2_obs=%.4f (p=%.4f)  B=%d (%d convergidas)\n",
            res_kuma$W2_obs, res_kuma$pW2, res_kuma$A2_obs, res_kuma$pA2,
            res_kuma$B, sum(res_kuma$conv_b)))
cat(sprintf("Beta:        W2_obs=%.4f (p=%.4f)  A2_obs=%.4f (p=%.4f)  B=%d (%d convergidas)\n",
            res_beta$W2_obs, res_beta$pW2, res_beta$A2_obs, res_beta$pA2,
            res_beta$B, sum(res_beta$conv_b)))
