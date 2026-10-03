## ARL: definicao operacional e controle.
## ARL0 = numero medio de observacoes ate o 1o alarme quando o processo esta SOB
## controle (queremos GRANDE: alarmes falsos raros).
## ARL1 = idem quando ha um deslocamento delta (queremos PEQUENO: deteccao rapida).
## Controla-se ARL0 escolhendo h. Aqui: k=0.5, residuos N(0,1) sob controle.
source("config_paths.R")  # Rodada 33: resolve dados/, resultados/, figuras/
set.seed(2026)
arl<-function(h,delta=0,k=0.5,B=6000,L=6000){
  rl<-integer(B)
  for(b in 1:B){C<-0;t<-0L
    repeat{t<-t+1L; C<-max(0, rnorm(1,mean=delta)-k+C)
      if(C>h){rl[b]<-t;break}; if(t>=L){rl[b]<-L;break}}}
  c(media=mean(rl), mediana=median(rl))
}
cat("=== ARL0: como h controla a taxa de alarme falso (delta=0) ===\n")
hs<-c(3.0,3.5,4.0,4.5,5.0)
t0<-sapply(hs,function(h) arl(h)["media"])
print(data.frame(h=hs, ARL0=round(t0,0), alarme_falso_a_cada=paste(round(t0/52,1),"anos")),row.names=FALSE)

cat("\n=== ARL1: poder de deteccao com h=4 (ARL0 ~ 340) ===\n")
ds<-c(0.5,1.0,1.5,2.0)
t1<-sapply(ds,function(d) arl(4.0,delta=d)["media"])
print(data.frame(deslocamento_delta=ds, ARL1=round(t1,1),
                 interpretacao=paste("detecta em ~",round(t1,1),"semanas")),row.names=FALSE)
saveRDS(list(hs=hs,arl0=t0,ds=ds,arl1=t1),"arl_result.rds")

## Figura: curva ARL0(h) e ARL1(delta)
cairo_pdf("Fig_arl.pdf", width=8.4, height=3.4, family="DejaVu Sans")
par(mfrow=c(1,2), mar=c(3.6,3.9,1.4,0.8), mgp=c(2.3,0.7,0), cex.axis=0.85, las=1)
plot(hs,t0,type="b",pch=16,xlab="Limite h",ylab=expression(ARL[0]),main="Alarme falso sob controle",cex.main=0.95)
abline(h=340,lty=3,col="grey50"); abline(v=4,lty=3,col="grey50")
plot(ds,t1,type="b",pch=16,xlab=expression(paste("Deslocamento ",delta)),ylab=expression(ARL[1]),main="Rapidez de detecção (h=4)",cex.main=0.95)
dev.off()
cat("Fig_arl gerada\n")
