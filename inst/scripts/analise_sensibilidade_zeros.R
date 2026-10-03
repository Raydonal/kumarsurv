###############################################################################
## ANALISE DE SENSIBILIDADE A SUBSTITUICAO DE ZEROS
## Varia epsilon (valor que substitui os zeros) e reajusta Kuma e beta para cada,
## medindo (i) ajuste Fase I (loglik, AIC) e (ii) calibracao dos limites Fase II.
## Objetivo: verificar se a VANTAGEM DE CALIBRACAO da Kuma e' robusta a epsilon,
## ou se e' artefato da escolha epsilon=1e-7. Responde ao parecer.
## Rodar: Rscript analise_sensibilidade_zeros.R
###############################################################################
source("copula_ml.R")
d <- read.csv("dda_platina_sem.csv", stringsAsFactors=FALSE)
y_raw <- d$incidence
X <- as.matrix(cbind(1, d$trend, d$cos26, d$sin26)); Z <- X
phaseI <- 300; I <- seq_len(phaseI); n <- nrow(d); II <- (phaseI+1):n
yII <- y_raw[II]                      # Fase II na escala ORIGINAL (zeros reais)

epsilons <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)
cat("=== SENSIBILIDADE A epsilon (substituicao de zeros) ===\n")
cat("Fase I: ajuste; Fase II: calibracao dos limites (alvo = tau)\n\n")
cat(sprintf("%-8s | %-22s | %-22s | %s\n","epsilon","Kuma loglik/AIC","Beta loglik/AIC","calib tau=0.95 (Kuma|beta)"))
cat(paste(rep("-",90),collapse=""),"\n")

fit_grade_tau <- function(y){
  # ajuste multi-tau validado, retorna limite calibrado em cada tau (Fase II)
  fits <- fit_copula_multi_tau(y[I], X[I,], Z[I,], p=1, q=1, family="kuma",
                               taus=c(0.50,0.75,0.90,0.95,0.99), verbose=FALSE)
  fits
}

res <- list()
for(eps in epsilons){
  y <- ifelse(y_raw==0, eps, y_raw)
  fk <- fit_copula(y[I], X[I,], Z[I,], p=1, q=1, family="kuma", tau=0.5)
  fb <- fit_copula(y[I], X[I,], Z[I,], p=1, q=1, family="beta", tau=0.5)
  # calibracao tau=0.95: Kuma estimada em 0.95 (validada) vs beta derivada da media
  fits <- fit_grade_tau(y)
  muk95 <- plogis(as.numeric(X[II,] %*% fits[["0.95"]]$par[1:ncol(X)]))
  mub  <- plogis(as.numeric(X[II,] %*% fb$par[1:ncol(X)]))
  phib <- exp(as.numeric(Z[II,] %*% fb$par[(ncol(X)+1):(ncol(X)+ncol(Z))]))
  qb95 <- qbeta(0.95, mub*phib, (1-mub)*phib)
  ck <- mean(yII <= muk95); cb <- mean(yII <= qb95)
  res[[as.character(eps)]] <- c(eps=eps, kloglik=fk$loglik, kaic=fk$aic,
                                bloglik=fb$loglik, baic=fb$aic, ck=ck, cb=cb)
  cat(sprintf("%-8.0e | %8.2f / %9.2f | %8.2f / %9.2f | %.3f | %.3f\n",
              eps, fk$loglik, fk$aic, fb$loglik, fb$aic, ck, cb))
}
cat("\nLeitura: se a coluna 'calib Kuma' permanece proxima de 0.95 e superior a\n")
cat("'calib beta' ao longo de todos os epsilon, a vantagem de calibracao NAO e'\n")
cat("artefato da escolha epsilon=1e-7 (o valor usado no artigo).\n")

M <- do.call(rbind, res)
cat(sprintf("\nAmplitude da calibracao Kuma tau=0.95 sobre epsilon: [%.3f, %.3f]\n",
            min(M[,"ck"]), max(M[,"ck"])))
cat(sprintf("Amplitude da calibracao beta tau=0.95 sobre epsilon:  [%.3f, %.3f]\n",
            min(M[,"cb"]), max(M[,"cb"])))
saveRDS(M, "sensibilidade_zeros.rds")
