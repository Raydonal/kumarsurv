## R47_valida_lib.R -- compara, variando so' o AR do ajuste adotado, a log-verossimilhanca da
## biblioteca corrigida (copula_ml.R), a da biblioteca anterior a Rodada 47
## (verificacao_R47/copula_ml_pre_R47.R, soma de 2000 pesos psi e Durbin-Levinson com L = 80)
## e a exata (Durbin-Levinson completo). Rodar a partir de scripts/.
source("copula_ml.R"); src_old <- new.env(); sys.source(file.path("verificacao_R47", "copula_ml_pre_R47.R"), envir = src_old)
R <- readRDS("reais_result_sem.rds"); I <- seq_len(R$phaseI); y <- R$y[I]; X <- R$X[I, ]; Z <- R$Z[I, ]; f <- R$fits[["kuma 1 1"]]
ll_new <- function(p) -nll_copula(p, y, X, Z, 1, 1, "kuma", 0.5); ll_old <- function(p) -src_old$nll_copula(p, y, X, Z, 1, 1, "kuma", 0.5)
ll_ex <- function(p) { mu <- plogis(as.numeric(X %*% p[1:4])); sh <- exp(as.numeric(Z %*% p[5:8])); eps <- qnorm(pmin(pmax(pkuma(y, mu, sh, .5), 1e-12), 1 - 1e-12))
  mom <- .one_step_moments(eps, .dl_recursion(as.numeric(ARMAacf(ar = p[9], ma = p[10], lag.max = length(y) - 1))))
  sum(dkuma(y, mu, sh, .5, log = TRUE)) - 0.5 * (sum(log(mom$s2)) + sum((eps - mom$m)^2 / mom$s2) - sum(eps^2)) }
for (a in c(0.5, 0.9, 0.948, 0.99, 0.996, 0.999, 0.9995)) { p <- f$par; p[9] <- a
  cat(sprintf("AR=%.4f  nova %.6f  antiga %.6f  exata %.6f\n", a, ll_new(p), ll_old(p), ll_ex(p))) }
t0 <- Sys.time(); for (i in 1:50) ll_new(f$par); t1 <- Sys.time(); for (i in 1:50) ll_old(f$par); t2 <- Sys.time()
cat(sprintf("tempo por avaliacao: nova %.1f ms, antiga %.1f ms\n", 1000 * as.numeric(t1 - t0) / 50, 1000 * as.numeric(t2 - t1) / 50))
## ponto de maior log-verossimilhanca encontrado pela multipartida de R47_perfil_ar1_global.R
## quando executado com a biblioteca ANTERIOR (AR = 0.9995): o valor alto era artefato do
## truncamento, e a verossimilhanca exata nesse ponto e' bem menor
p_esp <- c(-5.5035013447, -8.6834177186, -0.5870461821, -0.5196166085, -3.1459531466, -0.4549947368,
           0.0863077415, 0.0623517998, 0.9995346612, -0.8848084819)
cat(sprintf("\nponto espurio da biblioteca anterior (AR=%.4f): nova %.3f  antiga %.3f  exata %.3f\n",
            p_esp[9], ll_new(p_esp), ll_old(p_esp), ll_ex(p_esp)))
