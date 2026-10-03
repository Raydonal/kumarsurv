## R47_teste_kalman.R -- compara a forma quadratica pelo filtro de Kalman (com a variancia
## calculada, como antes da Rodada 47, por 2000 pesos psi) e o log-determinante de
## .logdet_arma com os exatos (Durbin-Levinson completo), na aplicacao e no desenho de Varin.
## Longe da fronteira (os dois casos) as vias coincidem.
source("copula_ml.R")
exato <- function(eps, ar, ma) {
  n <- length(eps); rho <- as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = n - 1))
  mom <- .one_step_moments(eps, .dl_recursion(rho))
  c(quad = sum((eps - mom$m)^2 / mom$s2), logdet = sum(log(mom$s2)))
}
kalman <- function(eps, ar, ma, ssinit = "Gardner1980") {
  n <- length(eps); psi <- ARMAtoMA(ar = ar, ma = ma, lag.max = 2000); g0 <- 1 + sum(psi^2)
  mod <- makeARIMA(phi = ar, theta = ma, Delta = numeric(0), SSinit = ssinit)
  kr <- KalmanRun(eps * sqrt(g0), mod, update = FALSE)
  c(quad = n * unname(kr$values[2]), logdet = .logdet_arma(ar, ma, n, 80))
}
## (a) aplicacao: ajuste Kuma(1,1) arquivado
R <- readRDS("reais_result_sem.rds"); f <- R$fits[["kuma 1 1"]]; X <- R$X; Z <- R$Z; kx <- ncol(X)
I <- seq_len(R$phaseI); y <- R$y[I]
mu <- plogis(as.numeric(X[I, ] %*% f$par[1:kx])); sh <- exp(as.numeric(Z[I, ] %*% f$par[(kx+1):(2*kx)]))
eps <- qnorm(pmin(pmax(pkuma(y, mu, sh, 0.5), 1e-12), 1 - 1e-12)); ar <- f$par[2*kx+1]; ma <- f$par[2*kx+2]
cat(sprintf("aplicacao: ar=%.4f ma=%.4f\n", ar, ma))
print(rbind(exato = exato(eps, ar, ma), kalman_Gardner = kalman(eps, ar, ma), kalman_Rossignol = kalman(eps, ar, ma, "Rossignol2011")), digits = 10)
## (b) desenho de Varin: uma replica simulada nos valores verdadeiros
set.seed(2024); N <- 52 * 9; time <- ((1:N) - 0.5 * N) / 100
Xv <- cbind(1, time, cos(2 * pi * (1:N) / 52), sin(2 * pi * (1:N) / 52))
e2 <- as.numeric(arima.sim(list(ar = c(1.5, -0.6), ma = -0.3), n = N)); e2 <- e2 / sd(e2)
cat("\nVarin: ar=(1.5,-0.6), ma=-0.3\n")
print(rbind(exato = exato(e2, c(1.5, -0.6), -0.3), kalman_Gardner = kalman(e2, c(1.5, -0.6), -0.3),
            kalman_Rossignol = kalman(e2, c(1.5, -0.6), -0.3, "Rossignol2011")), digits = 10)
## (c) curvatura em ar1 ao redor do ajuste da aplicacao: segunda diferenca numerica
g <- function(a, fun) { if (fun == "exato") v <- exato(eps, a, ma) else v <- kalman(eps, a, ma); 0.5 * (v["quad"] + v["logdet"]) }
h <- 1e-3
for (fun in c("exato", "kalman")) cat(sprintf("\nd2/dar1^2 da parte de copula (%s): %.2f", fun,
  (g(ar + h, fun) - 2 * g(ar, fun) + g(ar - h, fun)) / h^2))
cat("\n")
