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
REP <- 60; est <- matrix(NA, REP, K); set.seed(2025); ok <- 0
for(r in 1:REP){
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
print(out, row.names=FALSE)
cat("\nLeitura: vies pequeno em relacao ao DP indica recuperacao adequada dos\n")
cat("parametros no desenho efetivamente usado na aplicacao.\n")
