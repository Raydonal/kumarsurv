## simaplic_empates.R -- Rodada 47 (Quadro 4 do suplementar e Secao 5.2 do artigo). Por que, no desenho da aplicacao (Quadro 4,
## simulacao_desenho_aplicado_300.R), a estimativa do AR da copula fica junto a fronteira de
## estacionariedade, embora o processo gerador tenha AR = 0.948?
## O script de simulacao trunca as series em 1e-7 (pmax(y, 1e-7)). Com quantil mediano da
## ordem de 1e-4, cerca de um quarto dos valores simulados cai abaixo de 1e-7 e vira um empate
## em 1e-7, como os zeros substituidos por epsilon nos dados. Aqui as mesmas series (mesmas
## sementes, set.seed(2025 + r)) sao ajustadas com o mesmo procedimento, com e sem o truncamento
## (sem ele, os valores minusculos sao mantidos, limitados a 1e-300).
source("copula_ml.R")
zz <- file(file.path(DIR_RESULTADOS, "resultado_simaplic_empates.txt"), open = "wt"); sink(zz, split = TRUE)
R0 <- readRDS("reais_result_sem.rds"); tv <- R0$fits[["kuma 1 1"]]$par
N <- 300; s <- 1:N; time <- (s - mean(s)) / 100
X <- cbind(1, time, cos(2 * pi * s / 26), sin(2 * pi * s / 26)); Z <- X
fitfast <- function(y) {
  o0 <- optim(c(qlogis(median(y)), rep(0, 7)), function(p) nll_copula(p, y, X, Z, 0, 0, "kuma"), method = "BFGS", control = list(maxit = 300))
  optim(c(o0$par, 0.5, -0.3), function(p) nll_copula(p, y, X, Z, 1, 1, "kuma"), method = "BFGS", control = list(maxit = 400, reltol = 1e-9))
}
REP <- as.integer(Sys.getenv("REP", "50"))
res <- NULL
for (r in seq_len(REP)) {
  set.seed(2025 + r)
  y0 <- sim_copula(X, Z, beta = tv[1:4], gamma = tv[5:8], ar = tv[9], ma = tv[10], family = "kuma")
  oc <- try(fitfast(pmin(pmax(y0, 1e-7), 1 - 1e-7)), silent = TRUE)
  ou <- try(fitfast(pmin(pmax(y0, 1e-300), 1 - 1e-12)), silent = TRUE)
  res <- rbind(res, data.frame(rep = r, frac_empates = mean(y0 < 1e-7),
    ar_com_truncamento = if (inherits(oc, "try-error") || oc$convergence != 0) NA else oc$par[9],
    ar_sem_truncamento = if (inherits(ou, "try-error") || ou$convergence != 0) NA else ou$par[9]))
}
options(width = 160); print(res, digits = 4, row.names = FALSE)
cat(sprintf("\nfracao media de valores simulados abaixo de 1e-7 (empates apos o truncamento): %.3f\n", mean(res$frac_empates)))
for (v in c("ar_com_truncamento", "ar_sem_truncamento")) {
  a <- res[[v]]; cat(sprintf("%-20s convergidas %d de %d | mediana %.3f | fracao > 0.98: %.3f | fracao em [0.90, 0.98]: %.3f\n",
    v, sum(!is.na(a)), REP, median(a, na.rm = TRUE), mean(a > 0.98, na.rm = TRUE), mean(a >= 0.90 & a <= 0.98, na.rm = TRUE)))
}
sink(); close(zz)
saveRDS(res, file.path(DIR_RESULTADOS, "simaplic_empates_result.rds"))
