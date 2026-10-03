## R47_tab_varin.R -- Rodada 47: Quadro 8 (tab:mc-varin) com os erros-padrao
## recalculados por R47_ep_varin.R (passo de hessiana d = 1e-3) nas estimativas
## arquivadas de 14_mc_validacao_B300.R.
source("copula_ml.R")
m <- readRDS("mc_validacao_B300_result.rds")
## R47_ep_varin.R grava um arquivo por bloco de replicas (I0, I1): usa o de 1 a 300, se
## existir, ou os dois blocos de 150 com que a auditoria foi feita em paralelo
blocos <- if (file.exists(file.path(DIR_RESULTADOS, "R47_ep_varin_001_300.rds"))) "R47_ep_varin_001_300.rds" else
  c("R47_ep_varin_001_150.rds", "R47_ep_varin_151_300.rds")
L <- lapply(file.path(DIR_RESULTADOS, blocos), readRDS)
se_a <- do.call(rbind, lapply(L, `[[`, "se_arma")); se_i <- do.call(rbind, lapply(L, `[[`, "se_ind"))
stopifnot(nrow(se_a) == 300)
E <- m$est_arma; Ei <- m$est_ind
rob <- function(x) IQR(x, na.rm = TRUE) / 1.349
tab <- data.frame(param = m$tab$param, verdadeiro = m$tab$verdadeiro,
  media_arma = colMeans(E), dp_arma = apply(E, 2, sd), dprob_arma = apply(E, 2, rob),
  ep_arma_R47 = colMeans(se_a, na.rm = TRUE), ep_arma_mediana = apply(se_a, 2, median, na.rm = TRUE), ep_arma_publicado = m$tab$ep_arma,
  n_ep_ok = colSums(is.finite(se_a)),
  media_ind = c(colMeans(Ei), NA, NA, NA), dp_ind = c(apply(Ei, 2, sd), NA, NA, NA), dprob_ind = c(apply(Ei, 2, rob), NA, NA, NA),
  ep_ind_R47 = c(colMeans(se_i, na.rm = TRUE), NA, NA, NA), ep_ind_publicado = m$tab$ep_ind)
out <- file.path(DIR_RESULTADOS, "resultado_tab_varin_R47.txt"); zz <- file(out, "wt"); sink(zz, split = TRUE)
cat("=== Quadro 8 (desenho de Varin), EP recalculados com passo de hessiana d = 1e-3 ===\n")
options(width = 220); print(tab, digits = 3, row.names = FALSE)
atip <- rowSums(abs(sweep(E, 2, apply(E, 2, median))) > 5 * matrix(apply(E, 2, rob), nrow(E), ncol(E), byrow = TRUE)) > 0
cat(sprintf("\nreplicas com alguma estimativa a mais de 5 DP robustos da mediana: %d de %d\n", sum(atip), nrow(E)))
cat("razao DP/EP(R47):", sprintf("%.2f", tab$dp_arma / tab$ep_arma_R47), "\n")
cat("razao DProb/EP(R47):", sprintf("%.2f", tab$dprob_arma / tab$ep_arma_R47), "\n")
cat("sem as replicas atipicas, DP:", sprintf("%.3f", apply(E[!atip, ], 2, sd)), "\n")
sink(); close(zz)
saveRDS(list(tab = tab, se_arma = se_a, se_ind = se_i, atipicas = which(atip)), file.path(DIR_RESULTADOS, "tab_varin_R47.rds"))
