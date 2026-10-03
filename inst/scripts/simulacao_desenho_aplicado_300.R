###############################################################################
## SIMULACAO NO DESENHO DA APLICACAO
## ARMA(1,1) + harmonico SEMESTRAL (periodo 26), marginal Kumaraswamy
## Parametros verdadeiros = estimativas da Fase I (Secao 5), tornando o
## experimento diretamente informativo para a aplicacao.
## Rodar: Rscript simulacao_desenho_aplicado.R
###############################################################################
source("copula_ml.R")
R0 <- readRDS("reais_result_sem.rds"); f <- R0$fits[["kuma 1 1"]]
tv <- f$par                      # valores verdadeiros = estimativas da aplicacao
N  <- 300                        # mesmo tamanho da Fase I
s  <- 1:N; time <- (s-mean(s))/100
X  <- cbind(1, time, cos(2*pi*s/26), sin(2*pi*s/26)); Z <- X
colnames(X)<-colnames(Z)<-c("int","trend","cos26","sin26")
K  <- length(tv)
cat("Desenho: ARMA(1,1), harmonico semestral, n =",N,"\n")
cat("Valores verdadeiros (estimativas da Fase I):\n")
print(round(tv,3))

fitfast <- function(y){
  o0 <- optim(c(qlogis(median(y)),rep(0,7)), function(p) nll_copula(p,y,X,Z,0,0,"kuma"),
              method="BFGS", control=list(maxit=300))
  optim(c(o0$par, 0.5, -0.3), function(p) nll_copula(p,y,X,Z,1,1,"kuma"),
        method="BFGS", control=list(maxit=400, reltol=1e-9))
}
REP <- 300  # elevado de 60
## Rodada 41: a semeadura passa a ser POR REPLICA, set.seed(2025 + r), identica a
## de simulacao_checkpoint.R. Antes havia um unico set.seed(2025) e um fluxo
## aleatorio continuo, de modo que este script NAO reproduzia o Quadro 4 do
## suplementar (tab:simaplic), cujos numeros vieram da via incremental do
## checkpoint. Com a semeadura por replica as duas vias dao exatamente as mesmas
## 300 series, e a tabela passa a ser reproduzivel por uma unica execucao.
est <- matrix(NA, REP, K); ok <- 0
for(r in 1:REP){
  set.seed(2025 + r)
  y <- try(sim_copula(X, Z, beta=tv[1:4], gamma=tv[5:8], ar=tv[9], ma=tv[10],
                      family="kuma"), silent=TRUE)
  if(inherits(y,"try-error")) next
  y <- pmin(pmax(y,1e-7),1-1e-7)
  o <- try(fitfast(y), silent=TRUE)
  if(inherits(o,"try-error")||o$convergence!=0) next
  est[r,] <- o$par; ok <- ok+1
}
est <- est[complete.cases(est),,drop=FALSE]
cat(sprintf("\nReplicas convergidas: %d de %d\n\n", nrow(est), REP))
nm <- names(tv); if(is.null(nm)) nm <- paste0("p",1:K)
out <- data.frame(parametro=nm, verdadeiro=round(tv,3),
                  media=round(colMeans(est),3),
                  vies=round(colMeans(est)-tv,3),
                  vies_rel_pct=round(100*(colMeans(est)-tv)/abs(tv),1),
                  DP=round(apply(est,2,sd),3),
                  EQM=round(colMeans((est-matrix(tv,nrow(est),K,byrow=TRUE))^2),4))
## Rodada 47: desvio-padrao robusto (amplitude interquartil / 1.349) e replicas atipicas,
## como no Quadro 8; as estimativas por replica passam a ser gravadas.
out$DPr <- round(apply(est, 2, function(x) IQR(x) / 1.349), 3)
print(out, row.names=FALSE)
dpr <- apply(est, 2, function(x) IQR(x) / 1.349)
atip <- rowSums(abs(sweep(est, 2, apply(est, 2, median))) > 5 * matrix(dpr, nrow(est), K, byrow = TRUE)) > 0
cat(sprintf("\nreplicas com alguma estimativa a mais de 5 DPr da mediana: %d de %d\n", sum(atip), nrow(est)))
cat("DP sem essas replicas:", sprintf("%.3f", apply(est[!atip, , drop = FALSE], 2, sd)), "\n")
cat("media sem essas replicas:", sprintf("%.3f", colMeans(est[!atip, , drop = FALSE])), "\n")
cat(sprintf("AR(1): media %.3f, mediana %.3f, fracao de replicas com AR > 0.98: %.3f\n",
            mean(est[, K - 1]), median(est[, K - 1]), mean(est[, K - 1] > 0.98)))
saveRDS(list(est = est, tv = tv, tab = out, atipicas = which(atip)), "simulacao_aplicada_300_result.rds")
cat("\nLeitura: vies pequeno em relacao ao DP indica recuperacao adequada dos\n")
cat("parametros no desenho efetivamente usado na aplicacao.\n")
