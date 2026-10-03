###############################################################################
## ANALISE MULTI-QUANTIL — o argumento central do modelo quantilico
## Ajusta a marginal Kumaraswamy reparametrizada em varios niveis tau e mostra
## (i) como os coeficientes de regressao variam com tau (efeito heterogeneo)
## (ii) a calibracao do limite em cada tau (Fase II)
## Rodar: Rscript analise_multiquantil.R
###############################################################################
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds")
d<-R$d; X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI; n<-nrow(d)
I<-seq_len(phaseI); II<-(phaseI+1):n; yII<-d$incidence[II]
taus <- c(0.50,0.75,0.90,0.95,0.99)

cat("=================================================================\n")
cat(" (I) COEFICIENTES DO SUBMODELO DE LOCACAO POR NIVEL tau\n")
cat("     Kumaraswamy + copula ARMA(1,1), harmonico semestral, Fase I\n")
cat("=================================================================\n")
cat(sprintf("%-6s %10s %10s %10s %10s\n","tau","intercepto","tendencia","cos26","sin26"))
## Rodada 33: usar a selecao VALIDADA (monotonicidade dos quantis + busca
## multi-partida). fit_copula() puro converge para o otimo espurio em
## tau=0.90 e reproduzia a linha errada da tab:multitau.
fits_val <- fit_copula_multi_tau(y[I],X[I,],Z[I,],p=1,q=1,family="kuma",
                                 taus=taus, verbose=FALSE)
fits<-list()
for(tau0 in taus){
  f<-try(fits_val[[as.character(tau0)]],silent=TRUE)
  if(inherits(f,"try-error")){cat(sprintf("%-6.2f  falhou\n",tau0));next}
  fits[[as.character(tau0)]]<-f
  b<-f$par[1:ncol(X)]
  cat(sprintf("%-6.2f %10.3f %10.3f %10.3f %10.3f\n", tau0, b[1],b[2],b[3],b[4]))
}
cat("\nLeitura: se os coeficientes MUDAM com tau, o efeito das covariaveis nao e'\n")
cat("homogeneo ao longo da distribuicao — justamente o que um modelo de media\n")
cat("nao capta, e o argumento substantivo a favor da abordagem quantilica.\n")

cat("\n=================================================================\n")
cat(" (II) CALIBRACAO DO LIMITE EM CADA tau (Fase II; alvo = tau)\n")
cat("=================================================================\n")
cat(sprintf("%-6s %14s %14s\n","tau","Kuma@tau","desvio"))
for(tau0 in taus){
  f<-fits[[as.character(tau0)]]; if(is.null(f)) next
  mu<-plogis(as.numeric(X[II,]%*%f$par[1:ncol(X)]))
  cob<-mean(yII<=mu)
  cat(sprintf("%-6.2f %14.3f %+14.3f\n", tau0, cob, cob-tau0))
}
cat("\nLeitura: desvio proximo de zero indica limite bem calibrado no nivel pedido.\n")

cat("\n=================================================================\n")
cat(" (III) PERDA PINBALL NO PROPRIO NIVEL tau (Fase II) — menor e' melhor\n")
cat("  Kuma estimada EM tau  vs  Kuma estimada na mediana e extrapolada\n")
cat("=================================================================\n")
pinball<-function(q,tau0) mean(ifelse(yII>=q, tau0*(yII-q), (1-tau0)*(q-yII)))
f50<-fits[["0.5"]]
cat(sprintf("%-6s %14s %14s\n","tau","Kuma@tau","Kuma@0.50"))
for(tau0 in taus){
  f<-fits[[as.character(tau0)]]; if(is.null(f)||is.null(f50)) next
  q_at<-plogis(as.numeric(X[II,]%*%f$par[1:ncol(X)]))
  mu50<-plogis(as.numeric(X[II,]%*%f50$par[1:ncol(X)]))
  a50<-exp(as.numeric(Z[II,]%*%f50$par[(ncol(X)+1):(ncol(X)+ncol(Z))]))
  k50<-log(1-0.5)/log(1-mu50^a50)
  q_50<-(1-(1-tau0)^(1/k50))^(1/a50)
  cat(sprintf("%-6.2f %14.5f %14.5f\n", tau0, pinball(q_at,tau0), pinball(q_50,tau0)))
}

cat("\n=== (IV) Reajuste de tau=0.90 com valores iniciais do ajuste em 0.95 ===\n")
f95<-fits[["0.95"]]
if(!is.null(f95)){
  nll<-function(par) nll_copula(par,y=y[I],X=X[I,],Z=Z[I,],p=1,q=1,family="kuma",tau=0.90)
  o<-optim(f95$par, nll, method="BFGS", control=list(maxit=1500,reltol=1e-11))
  b<-o$par[1:ncol(X)]
  cat(sprintf("tau=0.90 reajustado: int=%.3f tend=%.3f cos=%.3f sin=%.3f (conv=%d)\n",
              b[1],b[2],b[3],b[4],o$convergence))
  mu<-plogis(as.numeric(X[II,]%*%b)); cob<-mean(yII<=mu)
  cat(sprintf("  calibracao: %.3f (alvo 0.90, desvio %+.3f)\n", cob, cob-0.90))
  cat(sprintf("  pinball: %.5f\n", mean(ifelse(yII>=mu, 0.90*(yII-mu), 0.10*(mu-yII)))))
}
