## English-label version of 07_fig_deteccao_SEM.R (final Fig_deteccao generator,
## overwrites 06's draft, per run_all.R ordering).
source("copula_ml.R")
if (!exists("DIR_FIGURAS_EN")) stop("DIR_FIGURAS_EN not set")
R<-readRDS("reais_result_sem.rds"); d<-R$d; phaseI<-R$phaseI; n<-nrow(d)
y_raw<-d$incidence
yr<-as.integer(substr(d$epiweek,1,4)); wk<-as.integer(substr(d$epiweek,7,8)); tnum<-yr+(wk-1)/52
ff<-readRDS("farrington_result.rds"); cs<-readRDS("cusum_result.rds")
af<-ff$alarms_idx; ac<-cs$flag
II<-(phaseI+1):n

cairo_pdf(file.path(DIR_FIGURAS_EN, "Fig_deteccao.pdf"), width=9, height=3.6, family="DejaVu Sans")
par(mar=c(3.4,3.9,0.6,1.8), mgp=c(2.3,0.7,0), cex.axis=0.9, las=1)
xr<-range(tnum[II])
plot(tnum[II], y_raw[II], type="h", col="grey70", lwd=1.4, xlim=xr, xaxs="i",
     ylim=c(-0.05,max(y_raw[II])*1.05), xlab="Year", ylab="ADD incidence", yaxt="n", xaxt="n")
axis(2, at=seq(0,0.20,0.05))
axis(1, at=2023:2026); axis(1, at=seq(2023.5, 2025.5, 1), labels=FALSE, tcl=-0.25)
rug_y1<--0.018; rug_y2<--0.038
segments(tnum[af], rug_y1-0.008, tnum[af], rug_y1+0.008, col="black", lwd=2.4)
segments(tnum[ac], rug_y2-0.008, tnum[ac], rug_y2+0.008, col="grey35", lwd=2.4)
text(xr[1], rug_y1, "Farrington", pos=2, cex=0.72, xpd=NA, offset=0.3)
text(xr[1], rug_y2, "CUSUM", pos=2, cex=0.72, xpd=NA, offset=0.3)
abline(h=0, col="grey85")
dev.off()
cat("Fig_deteccao (EN) ok\n")
