## English-label version of fig_quantis_corrigida.R (Figure S4 of the supplementary).
source("copula_ml.R")
if (!exists("DIR_FIGURAS_EN")) stop("DIR_FIGURAS_EN not set")
R<-readRDS("reais_result_sem.rds"); d<-R$d; X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI
I<-seq_len(phaseI)
yr<-as.integer(substr(d$epiweek,1,4)); wk<-as.integer(substr(d$epiweek,7,8))
y_raw<-d$incidence

fits_val<-fit_copula_multi_tau(y[I],X[I,],Z[I,],p=1,q=1,family="kuma",
                               taus=c(0.50,0.75,0.90,0.95,0.99), verbose=FALSE)
qf<-function(tt) plogis(as.numeric(X%*%fits_val[[as.character(tt)]]$par[1:ncol(X)]))
q50<-qf(0.5); q90<-qf(0.9); q95<-qf(0.95)
cat(sprintf("loglik tau=0.90: %.3f\n", fits_val[["0.9"]]$loglik))

sel<-which(yr==2022); eps<-1e-5
cairo_pdf(file.path(DIR_FIGURAS_EN, "Fig_quantis.pdf"), width=7.0, height=4.2, family="DejaVu Sans")
par(mar=c(3.6,5.0,1.0,1.0), mgp=c(2.6,0.7,0), cex.axis=0.9, las=1)
allv<-c(q50[sel],q90[sel],q95[sel],y_raw[sel]); allv<-allv[allv>0]
yl<-range(c(allv,eps))
plot(wk[sel], pmax(q95[sel],eps), type="l", lwd=1.8, log="y", ylim=yl,
     xlab="Epidemiological week", ylab="")
title(ylab="Marginal quantile (log scale)", line=3.7)
lines(wk[sel], pmax(q90[sel],eps), lwd=1.4, lty=2, col="#d95f02")
lines(wk[sel], pmax(q50[sel],eps), lwd=1.4, lty=3, col="#1b9e77")
obs<-y_raw[sel]; ok<-obs>0
points(wk[sel][ok], obs[ok], pch=1, cex=0.7, col="grey40")
legend("bottomright", bty="n", cex=0.8, ncol=2,
  legend=c(expression(tau==0.95), expression(tau==0.90), expression(tau==0.50), "observed incidence"),
  lty=c(1,2,3,NA), pch=c(NA,NA,NA,1), lwd=c(1.8,1.4,1.4,NA),
  col=c("black","#d95f02","#1b9e77","grey40"))
dev.off()
cat("Fig_quantis (EN) ok\n")
