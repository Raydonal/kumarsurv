source("copula_ml.R")
R<-readRDS("reais_result_sem.rds"); d<-R$d; phaseI<-R$phaseI; n<-nrow(d)
y_raw<-d$incidence
yr<-as.integer(substr(d$epiweek,1,4)); wk<-as.integer(substr(d$epiweek,7,8)); tnum<-yr+(wk-1)/52
ff<-readRDS("farrington_result.rds"); cs<-readRDS("cusum_result.rds")
af<-ff$alarms_idx; ac<-cs$flag
II<-(phaseI+1):n

cairo_pdf("Fig_deteccao.pdf", width=9, height=3.6, family="DejaVu Sans")
## Rodada 50: margem direita maior e eixo x em anos inteiros; o rotulo "2026.0",
## centrado na borda direita, saia cortado.
par(mar=c(3.4,3.9,0.6,1.8), mgp=c(2.3,0.7,0), cex.axis=0.9, las=1)
xr<-range(tnum[II])
plot(tnum[II], y_raw[II], type="h", col="grey70", lwd=1.4, xlim=xr, xaxs="i",
     ylim=c(-0.05,max(y_raw[II])*1.05), xlab="Ano", ylab="Incidência de DDA", yaxt="n", xaxt="n")
axis(2, at=seq(0,0.20,0.05))
axis(1, at=2023:2026); axis(1, at=seq(2023.5, 2025.5, 1), labels=FALSE, tcl=-0.25)
# faixas de alarme
rug_y1<--0.018; rug_y2<--0.038
segments(tnum[af], rug_y1-0.008, tnum[af], rug_y1+0.008, col="black", lwd=2.4)
segments(tnum[ac], rug_y2-0.008, tnum[ac], rug_y2+0.008, col="grey35", lwd=2.4)
text(xr[1], rug_y1, "Farrington", pos=2, cex=0.72, xpd=NA, offset=0.3)
text(xr[1], rug_y2, "CUSUM", pos=2, cex=0.72, xpd=NA, offset=0.3)
abline(h=0, col="grey85")
dev.off()
cat("Fig_deteccao redesenhada\n")
