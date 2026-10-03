## Validacao Monte Carlo da estimacao por verossimilhanca fechada (cópula
## gaussiana ARMA, marginal Kumaraswamy). Recuperacao de parametros no desenho
## ARMA(2,1) da Secao 4. Executar apos source("copula_ml.R").
source("copula_ml.R")
N<-52*9; time<-((1:N)-0.5*N)/100
X<-cbind(1,time,cos(2*pi*(1:N)/52),sin(2*pi*(1:N)/52)); Z<-X
colnames(X)<-colnames(Z)<-c("int","trend","cos","sin")
beta<-c(-2,0.5,-0.67,-0.22); gamma<-c(-2,0.1,-0.19,-0.06); ar<-c(1.5,-0.6); ma<--0.3
tv<-c(beta,gamma,ar,ma); K<-length(tv)
fitfast<-function(y){
  o0<-optim(c(qlogis(median(y)),rep(0,7)),function(p)nll_copula(p,y,X,Z,0,0,"kuma"),method="BFGS",control=list(maxit=200))
  optim(c(o0$par,rep(0.01,3)),function(p)nll_copula(p,y,X,Z,2,1,"kuma"),method="BFGS",control=list(maxit=200,reltol=1e-8))
}
R<-60; est<-matrix(NA,R,K); set.seed(2024)
for(i in 1:R){ y<-sim_copula(X,Z,beta,gamma,ar,ma,family="kuma")
  f<-try(suppressWarnings(fitfast(y)),silent=TRUE)
  if(!inherits(f,"try-error") && f$convergence==0) est[i,]<-f$par }
mn<-colMeans(est,na.rm=TRUE); sdv<-apply(est,2,sd,na.rm=TRUE)
nm<-c(paste0("mu.",colnames(X)),paste0("sh.",colnames(Z)),"ar1","ar2","ma1")
print(data.frame(param=nm,verdadeiro=round(tv,3),media=round(mn,3),
                 vies=round(mn-tv,3),SD=round(sdv,3)),row.names=FALSE)
