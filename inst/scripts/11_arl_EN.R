## English-label version of 11_arl.R, figure only (no saveRDS of arl_result.rds).
source("config_paths.R")
if (!exists("DIR_FIGURAS_EN")) stop("DIR_FIGURAS_EN not set")
set.seed(2026)
arl<-function(h,delta=0,k=0.5,B=6000,L=6000){
  rl<-integer(B)
  for(b in 1:B){C<-0;t<-0L
    repeat{t<-t+1L; C<-max(0, rnorm(1,mean=delta)-k+C)
      if(C>h){rl[b]<-t;break}; if(t>=L){rl[b]<-L;break}}}
  c(media=mean(rl), mediana=median(rl))
}
hs<-c(3.0,3.5,4.0,4.5,5.0)
t0<-sapply(hs,function(h) arl(h)["media"])
ds<-c(0.5,1.0,1.5,2.0)
t1<-sapply(ds,function(d) arl(4.0,delta=d)["media"])

cairo_pdf(file.path(DIR_FIGURAS_EN, "Fig_arl.pdf"), width=8.4, height=3.4, family="DejaVu Sans")
par(mfrow=c(1,2), mar=c(3.6,3.9,1.4,0.8), mgp=c(2.3,0.7,0), cex.axis=0.85, las=1)
plot(hs,t0,type="b",pch=16,xlab="Limit h",ylab=expression(ARL[0]),main="False alarm under control",cex.main=0.95)
abline(h=340,lty=3,col="grey50"); abline(v=4,lty=3,col="grey50")
plot(ds,t1,type="b",pch=16,xlab=expression(paste("Shift ",delta)),ylab=expression(ARL[1]),main="Detection speed (h=4)",cex.main=0.95)
dev.off()
cat("Fig_arl (EN) ok\n")
