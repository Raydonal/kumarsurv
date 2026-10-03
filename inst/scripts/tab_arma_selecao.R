source("copula_ml.R")
R<-readRDS("reais_result_sem.rds"); X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI; I<-seq_len(phaseI)
res<-list()
for(p in 0:3) for(q in 0:3){
  f<-tryCatch(fit_copula(y[I],X[I,],Z[I,],p=p,q=q,family="kuma",tau=0.5),error=function(e)NULL)
  if(is.null(f)) next
  r<-qresiduals(f)$r
  ac<-acf(r,plot=FALSE,lag.max=4)$acf[2:5]
  res[[paste(p,q)]]<-c(p=p,q=q,AIC=f$aic,ACF1=ac[1],ACF2=ac[2],ACF3=ac[3],ACF4=ac[4])
}
M<-do.call(rbind,res); M<-M[order(M[,"AIC"]),]
for(i in 1:nrow(M)) cat(sprintf("%d & %d & %.2f & %.2f & %.2f & %.2f & %.2f\n",
  M[i,"p"],M[i,"q"],M[i,"AIC"],M[i,"ACF1"],M[i,"ACF2"],M[i,"ACF3"],M[i,"ACF4"]))
