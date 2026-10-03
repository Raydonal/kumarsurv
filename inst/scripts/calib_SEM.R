## Calibracao dos limites de controle (Fase II) — harmonico SEMESTRAL
## IMPORTANTE: usa fit_copula_multi_tau() para obter os ajustes por nivel tau,
## que aplica o criterio de monotonicidade dos quantis + reinicializacao
## multi-partida. Isto e' necessario porque o ajuste padrao em tau=0.90 pode
## convergir para um otimo espurio (quantis acima do maximo observado,
## calibracao=1.000). O artigo reporta os valores destes ajustes validados.
source("copula_ml.R")
R<-readRDS("reais_result_sem.rds"); d<-R$d; X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI; n<-nrow(d)
I<-seq_len(phaseI); II<-(phaseI+1):n; yII<-d$incidence[II]
cat("=== CALIBRACAO DOS LIMITES DE CONTROLE (Fase II) — harmonico SEMESTRAL ===\n")
cat("Fracao de observacoes abaixo do limite de nivel tau (alvo = tau)\n\n")

# ajustes validados por nivel tau (monotonicidade + reinicializacao).
# Passa-se a grade COMPLETA (incluindo niveis baixos) porque a reinicializacao
# multi-partida precisa de ancoras bem-identificadas abaixo de 0.90 para escapar
# do otimo espurio nesse nivel; extraem-se depois os tres niveis de interesse.
fits <- fit_copula_multi_tau(y[I], X[I,], Z[I,], p=1, q=1, family="kuma",
                             taus=c(0.50,0.75,0.90,0.95,0.99), verbose=FALSE)
fb <- R$fits[["beta 1 1"]]
mub<-plogis(as.numeric(X[II,]%*%fb$par[1:ncol(X)]))
phib<-exp(as.numeric(Z[II,]%*%fb$par[(ncol(X)+1):(ncol(X)+ncol(Z))]))

cat(sprintf("%-6s %-14s %-14s\n","tau","Kuma@tau","Beta(media)"))
for(tau0 in c(0.90,0.95,0.99)){
  fk<-fits[[as.character(tau0)]]
  muk<-plogis(as.numeric(X[II,]%*%fk$par[1:ncol(X)]))
  qb<-qbeta(tau0, mub*phib, (1-mub)*phib)
  cat(sprintf("%-6.2f %-14.3f %-14.3f\n", tau0, mean(yII<=muk), mean(yII<=qb)))
}

## --- Identificacao fraca em tau=0.90: reportar AMBOS os otimos ---
## Rodada 41: os rotulos abaixo estavam TROCADOS e o comentario descrevia o
## estado anterior a Rodada 33. Os dois otimos NAO tem log-verossimilhanca
## equivalente: 2542.186 (reinic. de tau=0.75, calib 0.882) contra 2541.998
## (reinic. de tau=0.95, calib 0.941), diferenca de 0.188. O corpo da
## tab:multitau reporta o de MAIOR verossimilhanca (0.882), que e' o que
## fit_copula_multi_tau() seleciona pelo criterio de monotonicidade dos
## quantis; a nota de rodape registra o alternativo (0.941). Reproduzem-se
## aqui os dois para rastreabilidade completa.
cat("\n--- tau=0.90: identificacao fraca (dois otimos, loglik DIFERENTES) ---\n")
p95 <- fit_copula(y[I], X[I,], Z[I,], p=1, q=1, family="kuma", tau=0.95)$par
nll90 <- function(par) nll_copula(par, y=y[I], X=X[I,], Z=Z[I,], p=1, q=1, family="kuma", tau=0.90)
o_de95 <- optim(p95, nll90, method="BFGS", control=list(maxit=1500, reltol=1e-11))
mu90a <- plogis(as.numeric(X[II,]%*%o_de95$par[1:ncol(X)]))
cat(sprintf("  reinic. de tau=0.95: loglik=%.3f  calib=%.3f  (nota de rodape)\n",
            -o_de95$value, mean(yII<=mu90a)))
p75 <- fit_copula(y[I], X[I,], Z[I,], p=1, q=1, family="kuma", tau=0.75)$par
o_de75 <- optim(p75, nll90, method="BFGS", control=list(maxit=1500, reltol=1e-11))
mu90b <- plogis(as.numeric(X[II,]%*%o_de75$par[1:ncol(X)]))
cat(sprintf("  reinic. de tau=0.75: loglik=%.3f  calib=%.3f  (corpo da tabela; maior loglik)\n",
            -o_de75$value, mean(yII<=mu90b)))
