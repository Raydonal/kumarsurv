source("copula_ml.R")
R<-readRDS("reais_result_sem.rds"); d<-R$d; X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI; n<-nrow(d)
idx<-seq_len(phaseI); II<-(phaseI+1):n

# --- ajustar Kumaraswamy DIRETAMENTE no nivel tau (regressao quantilica) ---
# e comparar a previsao do tau-quantil um passo a frente na Fase II.
taus<-c(0.50,0.90,0.95,0.99)

pred_qforecast<-function(par,family,tau_fit){
  kx<-ncol(X); kz<-ncol(Z)
  mu<-plogis(as.numeric(X%*%par[1:kx])); sh<-exp(as.numeric(Z%*%par[(kx+1):(kx+kz)]))
  ar<-par[kx+kz+1]; ma<-par[kx+kz+2]
  if(family=="kuma"){Ft<-pkuma(y,mu,sh,tau_fit)} else {Ft<-pbetamp(y,mu,sh)}
  eps<-qnorm(pmin(pmax(Ft,1e-12),1-1e-12))
  mom<-.one_step_moments(eps,.dl_recursion(as.numeric(ARMAacf(ar=ar,ma=ma,lag.max=n-1))))
  list(mu=mu,sh=sh,m=mom$m,v=mom$s2)
}
# previsao do tau0-quantil de Y_t (um passo), para um modelo com marginal em tau_fit
qforecast<-function(pe,family,tau_fit,tau0){
  p<-pnorm(pe$m+sqrt(pe$v)*qnorm(tau0))
  if(family=="kuma") qkuma(p,pe$mu,pe$sh,tau_fit) else qbetamp(p,pe$mu,pe$sh)
}
pinball<-function(q,tau0,idx){ yy<-y[idx]; qq<-q[idx]; mean(ifelse(yy>=qq,tau0*(yy-qq),(1-tau0)*(qq-yy))) }

# beta (media) e kuma-mediana: ajustados uma vez
fb<-R$fits[["beta 1 1"]]; peb<-pred_qforecast(fb$par,"beta",0.5)
fk50<-R$fits[["kuma 1 1"]]; pek50<-pred_qforecast(fk50$par,"kuma",0.5)

cat("=== Perda pinball fora da amostra no nivel tau (Fase II) — menor e' melhor ===\n")
cat(sprintf("%-6s %12s %12s %12s\n","tau","Kuma@tau","Kuma@0.5","Beta(media)"))
## Rodada 41: usar fit_copula_multi_tau(), nao fit_copula() puro. Em tau=0.90 o
## ajuste puro cai no otimo espurio (armadilha 3 da secao 5 do CONTEXTO_RETOMADA)
## e este script devolvia 0.00682/0.00726 na linha tau=0.90, valor que contradiz
## a tab:multitau do artigo. A grade completa inclui 0.75 como ancora.
fits_val <- fit_copula_multi_tau(y[idx],X[idx,],Z[idx,],p=1,q=1,family="kuma",
                                 taus=c(0.50,0.75,0.90,0.95,0.99), verbose=FALSE)
tab<-data.frame()
for(tau0 in taus){
  # Kuma ajustada NO nivel tau0 (ajuste validado)
  fkt<-fits_val[[as.character(tau0)]]
  pekt<-pred_qforecast(fkt$par,"kuma",tau0)
  qk_at<-qforecast(pekt,"kuma",tau0,tau0)
  qk_50<-qforecast(pek50,"kuma",0.5,tau0)
  qb   <-qforecast(peb ,"beta",0.5,tau0)
  a<-pinball(qk_at,tau0,II); b<-pinball(qk_50,tau0,II); c<-pinball(qb,tau0,II)
  cat(sprintf("%-6.2f %12.5f %12.5f %12.5f\n",tau0,a,b,c))
  tab<-rbind(tab,data.frame(tau=tau0,Kuma_at=a,Kuma_med=b,Beta=c))
}
saveRDS(tab,"quantil_result.rds")
cat("\nVencedor por nivel:\n")
for(i in 1:nrow(tab)){ v<-c("Kuma@tau","Kuma@0.5","Beta")[which.min(c(tab$Kuma_at[i],tab$Kuma_med[i],tab$Beta[i]))]
  cat(sprintf("  tau=%.2f: %s\n",tab$tau[i],v)) }
