## Fig_quantis CORRIGIDA: (a) curvas dos ajustes de fit_copula_multi_tau() (desde a
## Rodada 50; ver a nota abaixo); (b) escala logaritmica no eixo y, pois os
## quantis cobrem varias ordens de grandeza (tau=0.5 ~ 1e-4 vs tau=0.95 ~ 1e-1).
source("copula_ml.R")
R<-readRDS("reais_result_sem.rds"); d<-R$d; X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI
I<-seq_len(phaseI)
yr<-as.integer(substr(d$epiweek,1,4)); wk<-as.integer(substr(d$epiweek,7,8))
y_raw<-d$incidence

## Rodada 50: as curvas vem de fit_copula_multi_tau(), como em 06_figuras_SEM.R e
## analise_multiquantil.R. O reinicio de tau = 0.90 a partir do ajuste em tau = 0.95,
## usado ate' o v13, leva ao SEGUNDO ponto estacionario (loglik 2541.998, intercepto
## -3.52), e nao ao reportado na Tabela 3 do artigo (2542.186, intercepto -4.29): a
## curva tau = 0.90 da figura ficava colada 'a de tau = 0.95 (razao media 0.86 em 2022,
## contra 0.42 no ajuste reportado).
fits_val<-fit_copula_multi_tau(y[I],X[I,],Z[I,],p=1,q=1,family="kuma",
                               taus=c(0.50,0.75,0.90,0.95,0.99), verbose=FALSE)
qf<-function(tt) plogis(as.numeric(X%*%fits_val[[as.character(tt)]]$par[1:ncol(X)]))
q50<-qf(0.5); q90<-qf(0.9); q95<-qf(0.95)
cat(sprintf("loglik tau=0.90: %.3f\n", fits_val[["0.9"]]$loglik))
cat(sprintf("ranges: q50 [%.2e,%.2e] q90 [%.2e,%.2e] q95 [%.2e,%.2e]\n",
   min(q50),max(q50),min(q90),max(q90),min(q95),max(q95)))

sel<-which(yr==2022); eps<-1e-5
cairo_pdf("Fig_quantis.pdf", width=7.0, height=4.2, family="DejaVu Sans")
## Rodada 50: margem esquerda maior e rotulo do eixo y afastado; os rotulos da escala
## log ("1e-02"), horizontais, sobrepunham-se ao titulo do eixo.
par(mar=c(3.6,5.0,1.0,1.0), mgp=c(2.6,0.7,0), cex.axis=0.9, las=1)
allv<-c(q50[sel],q90[sel],q95[sel],y_raw[sel]); allv<-allv[allv>0]
yl<-range(c(allv,eps))
plot(wk[sel], pmax(q95[sel],eps), type="l", lwd=1.8, log="y", ylim=yl,
     xlab="Semana epidemiológica", ylab="")
title(ylab="Quantil marginal (escala log)", line=3.7)
lines(wk[sel], pmax(q90[sel],eps), lwd=1.4, lty=2, col="#d95f02")
lines(wk[sel], pmax(q50[sel],eps), lwd=1.4, lty=3, col="#1b9e77")
obs<-y_raw[sel]; ok<-obs>0
points(wk[sel][ok], obs[ok], pch=1, cex=0.7, col="grey40")
legend("bottomright", bty="n", cex=0.8, ncol=2,
  legend=c(expression(tau==0.95), expression(tau==0.90), expression(tau==0.50), "incidência observada"),
  lty=c(1,2,3,NA), pch=c(NA,NA,NA,1), lwd=c(1.8,1.4,1.4,NA),
  col=c("black","#d95f02","#1b9e77","grey40"))
dev.off()
cat("Fig_quantis (log, ajuste multi-tau)\n")
