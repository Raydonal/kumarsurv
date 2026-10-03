## Exemplo simulado da banda de predicao: o que ela e', e ela cobre o nominal?
source("copula_ml.R")
set.seed(7)
n<-260; tt<-((1:n)-n/2)/100
X<-cbind(1,tt,cos(2*pi*(1:n)/52),sin(2*pi*(1:n)/52)); Z<-X
colnames(X)<-colnames(Z)<-c("int","trend","cos","sin")
beta<-c(qlogis(0.25),0.15,-0.5,-0.2); gamma<-c(log(2.5),0,0,0)
ar<-0.6; ma<-0.3
y<-sim_copula(X,Z,beta,gamma,ar,ma,family="kuma")
fit<-suppressWarnings(fit_copula(y,X,Z,p=1,q=1,family="kuma"))
par<-fit$par
mu<-plogis(as.numeric(X%*%par[1:4])); sh<-exp(as.numeric(Z%*%par[5:8]))
arh<-par[9]; mah<-par[10]
u<-pkuma(y,mu,sh); eps<-qnorm(pmin(pmax(u,1e-12),1-1e-12))
mom<-.one_step_moments(eps,.dl_recursion(as.numeric(ARMAacf(ar=arh,ma=mah,lag.max=n-1))))
m<-mom$m; s<-sqrt(mom$s2)
# banda preditiva de um passo, nivel 1-alpha: quantis gamma_lo/gamma_hi da preditiva
band<-function(g) sapply(1:n,function(i) qkuma(pnorm(m[i]+s[i]*qnorm(g)),mu[i],sh[i]))
lo<-band(0.025); hi<-band(0.975)
cat(sprintf("=== Banda preditiva de 95%% (um passo) em dados SIMULADOS ===\n"))
cat(sprintf("cobertura empirica: %.3f  (nominal 0.95)\n", mean(y>=lo & y<=hi)))
cat(sprintf("cobertura 80%%: %.3f (nominal 0.80)\n", {l<-band(0.10);h<-band(0.90); mean(y>=l & y<=h)}))
cat(sprintf("cobertura 50%%: %.3f (nominal 0.50)\n", {l<-band(0.25);h<-band(0.75); mean(y>=l & y<=h)}))
saveRDS(list(y=y,lo=lo,hi=hi,mu=mu,m=m,s=s,band=band),"banda_result.rds")

## Figura didatica: serie simulada + banda + mediana condicional
cairo_pdf("Fig_banda_sim.pdf", width=9, height=3.8, family="DejaVu Sans")
par(mar=c(3.6,3.9,0.8,1.0), mgp=c(2.3,0.7,0), cex.axis=0.9, las=1)
plot(1:n, y, type="n", xlab="Tempo (semanas)", ylab="Resposta simulada", ylim=c(0,max(hi)*1.02), xaxs="i")
polygon(c(1:n,rev(1:n)), c(hi,rev(lo)), col="grey88", border=NA)
lines(1:n, band(0.5), col="grey35", lwd=1.2)
points(1:n, y, pch=16, cex=0.45, col="black")
out<-which(y<lo|y>hi); points(out, y[out], pch=1, cex=1.1, col="black", lwd=1.3)
legend("topleft", bty="n", cex=0.8,
  legend=c("Banda preditiva de 95%","Mediana preditiva","Observação","Fora da banda"),
  fill=c("grey88",NA,NA,NA), border=c("grey88",NA,NA,NA),
  lty=c(NA,1,NA,NA), pch=c(NA,NA,16,1), col=c(NA,"grey35","black","black"))
dev.off()
cat("Fig_banda_sim gerada\n")
