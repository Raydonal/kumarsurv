source("copula_ml.R")
R <- readRDS("reais_result_sem.rds")
d<-R$d; X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI; n<-nrow(d)
y_raw<-d$incidence
fk<-R$fits[["kuma 1 1"]]; fb<-R$fits[["beta 1 1"]]

pred_eval<-function(fit,family){
  par<-fit$par; kx<-ncol(X); kz<-ncol(Z)
  mu<-plogis(as.numeric(X%*%par[1:kx])); sh<-exp(as.numeric(Z%*%par[(kx+1):(kx+kz)]))
  ar<-par[kx+kz+1]; ma<-par[kx+kz+2]
  if(family=="kuma"){ Ft<-pkuma(y,mu,sh); ft<-dkuma(y,mu,sh); Qf<-function(u,i) qkuma(u,mu[i],sh[i]) }
  else { Ft<-pbetamp(y,mu,sh); ft<-dbetamp(y,mu,sh); Qf<-function(u,i) qbetamp(u,mu[i],sh[i]) }
  eps<-qnorm(pmin(pmax(Ft,1e-12),1-1e-12))
  rho<-as.numeric(ARMAacf(ar=ar,ma=ma,lag.max=n-1)); mom<-.one_step_moments(eps,.dl_recursion(rho))
  m<-mom$m; v<-mom$s2; z<-(eps-m)/sqrt(v)
  logscore <- -(dnorm(z,log=TRUE)-0.5*log(v)+log(ft)-dnorm(eps,log=TRUE))
  predq<-function(tau) sapply(seq_len(n),function(i) Qf(pnorm(m[i]+sqrt(v[i])*qnorm(tau)),i))
  list(logscore=logscore, predq=predq, m=m, v=v)
}
ek<-pred_eval(fk,"kuma"); eb<-pred_eval(fb,"beta")

II<-(phaseI+1):n                         # Fase II
pos<-II[y_raw[II]>0]                      # semanas positivas na Fase II
# rotulo de surto: semanas sinalizadas pelo Farrington (referencia)
ff<-readRDS("farrington_result.rds"); outbreak<-ff$alarms_idx
nonout<-setdiff(II, outbreak)

pinball<-function(e,tau,idx){ q<-e$predq(tau)[idx]; yy<-y[idx]; mean(ifelse(yy>=q,tau*(yy-q),(1-tau)*(q-yy))) }
taus<-c(0.90,0.95,0.975,0.99)
cat("=== Escore log preditivo fora da amostra (Fase II) — menor e' melhor ===\n")
cat(sprintf("  todas as semanas:   Kuma=%.3f  Beta=%.3f\n", mean(ek$logscore[II]), mean(eb$logscore[II])))
cat(sprintf("  semanas positivas:  Kuma=%.3f  Beta=%.3f\n", mean(ek$logscore[pos]), mean(eb$logscore[pos])))
cat("\n=== Perda pinball na cauda superior (Fase II, todas) — menor e' melhor ===\n")
for(tau in taus) cat(sprintf("  tau=%.3f:  Kuma=%.5f  Beta=%.5f  (%s)\n",tau,
   pinball(ek,tau,II),pinball(eb,tau,II), ifelse(pinball(ek,tau,II)<pinball(eb,tau,II),"Kuma","Beta")))
# CRPS aproximado por integral da pinball
grid<-seq(0.01,0.99,0.01)
crps_k<-2*mean(sapply(grid,function(t)pinball(ek,t,II))); crps_b<-2*mean(sapply(grid,function(t)pinball(eb,t,II)))
cat(sprintf("\n=== CRPS medio fora da amostra (Fase II) — menor e' melhor ===\n  Kuma=%.5f  Beta=%.5f  (%s)\n",
   crps_k,crps_b,ifelse(crps_k<crps_b,"Kuma","Beta")))
cat("\n=== Cobertura da banda superior nas semanas SEM surto (Fase II) ===\n")
for(tau in c(0.95,0.99)){ qk<-ek$predq(tau)[nonout]; qb<-eb$predq(tau)[nonout]
  cat(sprintf("  nivel %.0f%%:  Kuma cobre %.1f%%  Beta cobre %.1f%%  (alvo %.0f%%)\n",
    100*tau, 100*mean(y[nonout]<=qk), 100*mean(y[nonout]<=qb), 100*tau)) }
