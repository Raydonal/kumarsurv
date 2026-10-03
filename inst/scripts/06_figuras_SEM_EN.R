## English-label version of 06_figuras_SEM.R, Fig_carta panel only (Figure 3 of the
## paper). No saveRDS: this script reproduces the plot inputs, it does not
## replace the canonical resultados/calib_result.rds.
source("copula_ml.R")
if (!exists("DIR_FIGURAS_EN")) stop("DIR_FIGURAS_EN not set")
R<-readRDS("reais_result_sem.rds"); d<-R$d; X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI; n<-nrow(d)
y_raw<-d$incidence; idx<-seq_len(phaseI); II<-(phaseI+1):n
yr<-as.integer(substr(d$epiweek,1,4)); wk<-as.integer(substr(d$epiweek,7,8)); tnum<-yr+(wk-1)/52

## Same multi-start grid as 06_figuras_SEM.R: fit_copula_multi_tau(), never
## bare fit_copula(), so tau=0.90 lands on the reported stationary point.
taus_grade <- c(0.50,0.75,0.90,0.95,0.99)
fits_val <- fit_copula_multi_tau(y[idx],X[idx,],Z[idx,],p=1,q=1,family="kuma",
                                 taus=taus_grade, verbose=FALSE)
taus<-c(0.5,0.9,0.95); qfit<-list()
for(tt in taus){ f<-fits_val[[as.character(tt)]]
  qfit[[as.character(tt)]]<-plogis(as.numeric(X%*%f$par[1:4])) }
q50<-qfit[["0.5"]]; q90<-qfit[["0.9"]]; q95<-qfit[["0.95"]]
ncruz <- sum(q90<q50)+sum(q95<q90)
cat("crossings:", ncruz, "\n")
if(ncruz > 0) warning("non-monotone conditional quantiles in tau: ", ncruz,
                      " crossings, check the fit before using the figure")

cairo_pdf(file.path(DIR_FIGURAS_EN, "Fig_carta.pdf"), width=9, height=4.2, family="DejaVu Sans")
par(mar=c(3.6,3.8,1.2,1.0), mgp=c(2.3,0.7,0), cex.axis=0.9, las=1)
plot(tnum, y_raw, type="h", col="grey72", lwd=1, xlab="Year", ylab="ADD incidence",
     ylim=c(0,max(y_raw)*1.02), xaxs="i")
polygon(c(tnum,rev(tnum)), c(q95,rep(0,n)), col=rgb(0,0,0,0.05), border=NA)
abline(v=tnum[phaseI], lty=3, col="grey40")
lines(tnum, q50, col="grey45", lwd=1.2)
lines(tnum, q90, col="grey20", lwd=1.1, lty=2)
lines(tnum, q95, col="black", lwd=1.6)
exc<-II[y_raw[II]>q95[II]]; points(tnum[exc], y_raw[exc], pch=19, cex=0.65)
legend("topleft", bty="n", cex=0.82,
  legend=c("Weekly incidence","Marginal median (τ=0.5)","Control limit τ=0.90","Control limit τ=0.95","Exceeds the limit"),
  lty=c(1,1,2,1,NA), pch=c(NA,NA,NA,NA,19), lwd=c(1,1.2,1.1,1.6,NA),
  col=c("grey72","grey45","grey20","black","black"))
mtext("Phase I ", side=3, line=0.1, at=tnum[phaseI], adj=1, cex=0.8, col="grey40")
mtext(" Phase II", side=3, line=0.1, at=tnum[phaseI], adj=0, cex=0.8, col="grey40")
dev.off()
cat("Fig_carta (EN) ok\n")
