source("copula_ml.R")
n<-250; set.seed(7)
tt<-((1:n)-0.5*n)/100
X<-cbind(1,tt,cos(2*pi*(1:n)/52),sin(2*pi*(1:n)/52)); Z<-cbind(1); 
Xc<-X; Zc<-matrix(1,n,1)
beta<-c(qlogis(0.45),0.3,-0.5,-0.2)   # mediana/media espalhada em (0,1)
ar<-0.5; ma<-0.3
crps<-function(fit,family,y){
  par<-fit$par; kx<-ncol(Xc); kz<-1
  mu<-plogis(as.numeric(Xc%*%par[1:kx])); sh<-exp(as.numeric(Zc%*%par[(kx+1):(kx+kz)]))
  arh<-par[kx+kz+1]; mah<-par[kx+kz+2]
  if(family=="kuma"){Ft<-pkuma(y,mu,sh);Qf<-function(u,i)qkuma(u,mu[i],sh[i])} else {Ft<-pbetamp(y,mu,sh);Qf<-function(u,i)qbetamp(u,mu[i],sh[i])}
  eps<-qnorm(pmin(pmax(Ft,1e-12),1-1e-12))
  rho<-as.numeric(ARMAacf(ar=arh,ma=mah,lag.max=n-1)); mom<-.one_step_moments(eps,.dl_recursion(rho))
  m<-mom$m; v<-mom$s2; grid<-seq(0.05,0.95,0.05)
  mean(sapply(grid,function(tau){q<-sapply(1:n,function(i)Qf(pnorm(m[i]+sqrt(v[i])*qnorm(tau)),i)); mean(ifelse(y>=q,tau*(y-q),(1-tau)*(q-y)))}))*2
}
fitfast<-function(y,fam){
  o0<-optim(c(qlogis(median(y)),rep(0,ncol(Xc)-1+ncol(Zc))),function(p)nll_copula(p,y,Xc,Zc,0,0,fam),method="BFGS",control=list(maxit=150))
  optim(c(o0$par,0.01,0.01),function(p)nll_copula(p,y,Xc,Zc,1,1,fam),method="BFGS",control=list(maxit=150,reltol=1e-8))
}
R<-50
res<-array(NA,c(R,2,2),dimnames=list(NULL,c("DGP_kuma","DGP_beta"),c("fit_kuma","fit_beta")))
gam_k<-log(2.0)   # forma Kuma alpha=2 (concava, assimetrica)
gam_b<-log(12)    # precisao beta
for(i in 1:R){
  yk<-sim_copula(Xc,Zc,beta,gam_k,ar,ma,family="kuma")
  yb<-sim_copula(Xc,Zc,beta,gam_b,ar,ma,family="beta")
  for(dat in list(list(y=yk,d="DGP_kuma"),list(y=yb,d="DGP_beta"))){
    fk<-try(suppressWarnings(fitfast(dat$y,"kuma")),silent=TRUE)
    fb<-try(suppressWarnings(fitfast(dat$y,"beta")),silent=TRUE)
    if(!inherits(fk,"try-error")&&fk$convergence==0) res[i,dat$d,"fit_kuma"]<-crps(fk,"kuma",dat$y)
    if(!inherits(fb,"try-error")&&fb$convergence==0) res[i,dat$d,"fit_beta"]<-crps(fb,"beta",dat$y)
  }
}
cat("=== CRPS medio (menor e' melhor); linhas=verdade, colunas=modelo ajustado ===\n")
M<-apply(res,c(2,3),mean,na.rm=TRUE); print(round(M,5))
cat("\nGanho relativo do modelo CORRETO vs o incorreto:\n")
cat(sprintf("  DGP Kumaraswamy: beta perde %.1f%% em CRPS\n",100*(M["DGP_kuma","fit_beta"]/M["DGP_kuma","fit_kuma"]-1)))
cat(sprintf("  DGP Beta:        kuma perde %.1f%% em CRPS\n",100*(M["DGP_beta","fit_kuma"]/M["DGP_beta","fit_beta"]-1)))
saveRDS(M,"misspec_result.rds")
