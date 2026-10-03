###############################################################################
## Simulacao no desenho da aplicacao, com CHECKPOINT em disco.
## Cada chamada roda N_BATCH novas replicas e acrescenta ao checkpoint salvo,
## permitindo retomar entre chamadas (robusto a reinicios do ambiente que
## matam processos em segundo plano).
## Uso: Rscript simulacao_checkpoint.R <N_BATCH>
###############################################################################
source("copula_ml.R")
args <- commandArgs(trailingOnly = TRUE)
N_BATCH <- if (length(args) >= 1) as.integer(args[1]) else 30
CKPT <- "sim_checkpoint.rds"

R0 <- readRDS("reais_result_sem.rds"); f <- R0$fits[["kuma 1 1"]]
tv <- f$par
N  <- 300; s <- 1:N; time <- (s-mean(s))/100
X  <- cbind(1, time, cos(2*pi*s/26), sin(2*pi*s/26)); Z <- X
colnames(X)<-colnames(Z)<-c("int","trend","cos26","sin26")
K  <- length(tv)

if (file.exists(CKPT)) {
  ck <- readRDS(CKPT)
  est <- ck$est; done <- ck$done; seed_ctr <- ck$seed_ctr
  cat(sprintf("Retomando checkpoint: %d replicas ja concluidas.\n", done))
} else {
  est <- matrix(NA, 0, K); done <- 0; seed_ctr <- 0
  cat("Iniciando checkpoint novo.\n")
}

fitfast <- function(y){
  o0 <- optim(c(qlogis(median(y)),rep(0,7)), function(p) nll_copula(p,y,X,Z,0,0,"kuma"),
              method="BFGS", control=list(maxit=300))
  optim(c(o0$par, 0.5, -0.3), function(p) nll_copula(p,y,X,Z,1,1,"kuma"),
        method="BFGS", control=list(maxit=400, reltol=1e-9))
}

novas <- matrix(NA, N_BATCH, K); ok <- 0
for (b in 1:N_BATCH) {
  seed_ctr <- seed_ctr + 1
  set.seed(2025 + seed_ctr)
  y <- try(sim_copula(X, Z, beta=tv[1:4], gamma=tv[5:8], ar=tv[9], ma=tv[10], family="kuma"),
           silent = TRUE)
  if (inherits(y, "try-error")) next
  y <- pmin(pmax(y,1e-7),1-1e-7)
  o <- try(fitfast(y), silent = TRUE)
  if (inherits(o, "try-error") || o$convergence != 0) next
  ok <- ok + 1; novas[ok, ] <- o$par
}
novas <- novas[seq_len(ok), , drop = FALSE]
est <- rbind(est, novas); done <- done + N_BATCH
saveRDS(list(est = est, done = done, seed_ctr = seed_ctr), CKPT)

cat(sprintf("Bloco concluido: %d/%d replicas tentadas neste bloco convergiram.\n", ok, N_BATCH))
cat(sprintf("Total acumulado: %d tentativas, %d convergidas.\n", done, nrow(est)))

if (nrow(est) >= 30) {
  nm <- names(tv); if (is.null(nm)) nm <- paste0("p", 1:K)
  out <- data.frame(parametro = nm, verdadeiro = round(tv, 3),
                    media = round(colMeans(est), 3),
                    vies = round(colMeans(est) - tv, 3),
                    DP = round(apply(est, 2, sd), 3),
                    EQM = round(colMeans((est - matrix(tv, nrow(est), K, byrow = TRUE))^2), 4))
  cat("\n--- resumo parcial ate agora ---\n")
  print(out, row.names = FALSE)
}
