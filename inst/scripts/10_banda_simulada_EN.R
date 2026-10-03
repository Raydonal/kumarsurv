## English-label version of 10_banda_simulada.R, figure only (no saveRDS of
## banda_result.rds).
source("copula_ml.R")
if (!exists("DIR_FIGURAS_EN")) stop("DIR_FIGURAS_EN not set")
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
band<-function(g) sapply(1:n,function(i) qkuma(pnorm(m[i]+s[i]*qnorm(g)),mu[i],sh[i]))
lo<-band(0.025); hi<-band(0.975)

cairo_pdf(file.path(DIR_FIGURAS_EN, "Fig_banda_sim.pdf"), width=9, height=3.8, family="DejaVu Sans")
par(mar=c(3.6,3.9,0.8,1.0), mgp=c(2.3,0.7,0), cex.axis=0.9, las=1)
plot(1:n, y, type="n", xlab="Time (weeks)", ylab="Simulated response", ylim=c(0,max(hi)*1.02), xaxs="i")
polygon(c(1:n,rev(1:n)), c(hi,rev(lo)), col="grey88", border=NA)
lines(1:n, band(0.5), col="grey35", lwd=1.2)
points(1:n, y, pch=16, cex=0.45, col="black")
out<-which(y<lo|y>hi); points(out, y[out], pch=1, cex=1.1, col="black", lwd=1.3)
legend("topleft", bty="n", cex=0.8,
  legend=c("95% predictive band","Predictive median","Observation","Outside the band"),
  fill=c("grey88",NA,NA,NA), border=c("grey88",NA,NA,NA),
  lty=c(NA,1,NA,NA), pch=c(NA,NA,16,1), col=c(NA,"grey35","black","black"))
dev.off()
cat("Fig_banda_sim (EN) ok\n")
