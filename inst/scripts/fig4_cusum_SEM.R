## Regenera a Figura da carta CUSUM (Fig4) com o modelo SEMESTRAL
source("copula_ml.R")
R<-readRDS("reais_result_sem.rds"); d<-R$d; phaseI<-R$phaseI; n<-nrow(d)
C<-readRDS("cusum_result.rds")
yr<-as.integer(substr(d$epiweek,1,4)); wk<-as.integer(substr(d$epiweek,7,8))
tnum<-yr+(wk-1)/52
h<-4
cairo_pdf("Fig4.pdf", width=9, height=5.2, family="DejaVu Sans")
par(mfrow=c(2,1), mar=c(3.2,4.2,1.3,1.0), mgp=c(2.4,0.7,0), las=1, cex.axis=0.9)
## painel superior: incidencia + separacao Fase I / II
plot(tnum, d$incidence, type="h", col="grey45", lwd=1.1, xlab="", ylab="Incidência semanal", xaxs="i")
abline(v=tnum[phaseI], col="steelblue", lwd=2, lty=2)
## Rodada 50: rotulos na margem superior; dentro do grafico, "Fase II" colidia com
## a barra do pico do surto.
mtext(" Fase II", side=3, line=0.05, at=tnum[phaseI], adj=0, col="steelblue", cex=0.9)
mtext("Fase I ", side=3, line=0.05, at=tnum[phaseI], adj=1, col="steelblue", cex=0.9)
## painel inferior: estatistica CUSUM
cp<-C$cp; flag<-C$flag
plot(tnum[seq_along(cp)], cp, type="l", col="grey25", lwd=1.2,
     xlab="Ano", ylab=expression(C[t]^"+"), xaxs="i", ylim=c(0,max(cp,h)*1.05))
abline(h=h, col="firebrick", lwd=2, lty=2)
text(min(tnum) + 0.30, h, paste0("h = ",h), col="firebrick", pos=3, cex=0.9)
if(length(flag)) points(tnum[flag], cp[flag], pch=19, col="firebrick", cex=0.8)
dev.off()
cat("Fig4 (CUSUM) regenerada com modelo semestral\n")
