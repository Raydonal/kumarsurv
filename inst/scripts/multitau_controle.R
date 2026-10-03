## multitau_controle.R -- Rodada 47: colunas da tab:multitau avaliadas sob controle.
## Mesmos ajustes de analise_multiquantil.R (fit_copula_multi_tau); acrescenta a
## cobertura na Fase I e na Fase II sob controle e a perda pinball restrita as
## semanas sob controle, ao lado da perda em todas as semanas da Fase II.
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); d <- R$d; X <- R$X; Z <- R$Z; y <- R$y; phaseI <- R$phaseI; n <- nrow(d)
I <- seq_len(phaseI); II <- (phaseI + 1):n; yr <- d$incidence
F_out <- readRDS("farrington_result.rds")$alarms_idx; C_out <- readRDS("cusum_result.rds")$flag
ctrl <- setdiff(II, union(F_out, C_out)); taus <- c(0.50, 0.75, 0.90, 0.95, 0.99)
fits <- fit_copula_multi_tau(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma", taus = taus, verbose = FALSE)
f50 <- fits[["0.5"]]; kx <- ncol(X)
mu50 <- plogis(as.numeric(X %*% f50$par[1:kx])); sh50 <- exp(as.numeric(Z %*% f50$par[(kx + 1):(2 * kx)]))
pin <- function(q, t0, ix) mean(ifelse(yr[ix] >= q[ix], t0 * (yr[ix] - q[ix]), (1 - t0) * (q[ix] - yr[ix])))
out <- file.path(DIR_RESULTADOS, "resultado_multitau_controle.txt"); zz <- file(out, "wt"); sink(zz, split = TRUE)
cat("=== tab:multitau avaliada sob controle (Fase II sem as 25 semanas do surto) ===\n")
cat(sprintf("%-5s | %-7s %-7s %-7s | %-9s %-9s | %-9s %-9s\n", "tau", "cobI", "cobCtrl", "cobII", "pinII@t", "pinII@.5", "pinCt@t", "pinCt@.5"))
for (t0 in taus) {
  q_dir <- plogis(as.numeric(X %*% fits[[as.character(t0)]]$par[1:kx]))
  q_med <- qkuma(rep(t0, n), mu50, sh50, 0.5)
  cat(sprintf("%-5.2f | %.3f   %.3f   %.3f   | %.5f   %.5f   | %.5f   %.5f\n", t0,
              mean(yr[I] <= q_dir[I]), mean(yr[ctrl] <= q_dir[ctrl]), mean(yr[II] <= q_dir[II]),
              pin(q_dir, t0, II), pin(q_med, t0, II), pin(q_dir, t0, ctrl), pin(q_med, t0, ctrl)))
}
sink(); close(zz)
