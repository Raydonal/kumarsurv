## Tabela de simulacao comparando Kumaraswamy e beta:
## recuperacao (vies/EQM) do estimador sob cada marginal, mesmo desenho.
source("copula_ml.R")
n<-250; tt<-((1:n)-n/2)/100
X<-cbind(1,tt,cos(2*pi*(1:n)/52),sin(2*pi*(1:n)/52)); Z<-matrix(1,n,1)
colnames(X)<-c("int","trend","cos","sin"); colnames(Z)<-"int"
beta<-c(qlogis(0.30),0.25,-0.45,-0.18)
gk<-log(2.0); gb<-log(15)
ar<-0.6; ma<-0.3
fitq<-function(y,fam){
  o0<-optim(c(qlogis(median(y)),rep(0,4)),function(p)nll_copula(p,y,X,Z,0,0,fam),method="BFGS",control=list(maxit=300))
  optim(c(o0$par,0.01,0.01),function(p)nll_copula(p,y,X,Z,1,1,fam),method="BFGS",control=list(maxit=400,reltol=1e-9))
}
R<-120
run<-function(fam,gtrue){
  tv<-c(beta,gtrue,ar,ma); K<-length(tv); est<-matrix(NA,R,K)
  set.seed(99)
  for(i in 1:R){ y<-sim_copula(X,Z,beta,gtrue,ar,ma,family=fam)
    f<-try(suppressWarnings(fitq(y,fam)),silent=TRUE)
    if(!inherits(f,"try-error")&&f$convergence==0) est[i,]<-f$par }
  mn<-colMeans(est,na.rm=TRUE); sd<-apply(est,2,sd,na.rm=TRUE)
  data.frame(par=c("mu:int","mu:trend","mu:cos","mu:sin","forma","ar1","ma1"),
             verdadeiro=round(tv,3), media=round(mn,3),
             vies=round(mn-tv,3), EQM=round((mn-tv)^2+sd^2,4),
             conv=sum(!is.na(est[,1])))
}
cat("=== Kumaraswamy (marginal correta) ===\n"); rk<-run("kuma",gk); print(rk,row.names=FALSE)
cat("\n=== Beta (marginal correta) ===\n");        rb<-run("beta",gb); print(rb,row.names=FALSE)
saveRDS(list(kuma=rk,beta=rb),"simcomp_result.rds")
