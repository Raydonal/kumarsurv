source("copula_ml.R")
R<-readRDS("reais_result_sem.rds"); d<-R$d; X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI; n<-nrow(d)
y_raw<-d$incidence; idx<-seq_len(phaseI); II<-(phaseI+1):n
yr<-as.integer(substr(d$epiweek,1,4)); wk<-as.integer(substr(d$epiweek,7,8)); tnum<-yr+(wk-1)/52

## Rodada 41: usar fit_copula_multi_tau(), nao fit_copula() puro. Em tau=0.90
## o ajuste puro converge para o otimo espurio (calibracao 1.000, quantis acima
## do maximo observado), o que punha a curva tau=0.90 ACIMA da tau=0.95 em todos
## os 470 instantes e a deixava fora da area do grafico na Fig_carta. A grade
## completa e' necessaria: a busca multi-partida precisa das ancoras abaixo de
## 0.90 para escapar daquela bacia.
taus_grade <- c(0.50,0.75,0.90,0.95,0.99)
fits_val <- fit_copula_multi_tau(y[idx],X[idx,],Z[idx,],p=1,q=1,family="kuma",
                                 taus=taus_grade, verbose=FALSE)
taus<-c(0.5,0.9,0.95); qfit<-list()
for(tt in taus){ f<-fits_val[[as.character(tt)]]
  qfit[[as.character(tt)]]<-plogis(as.numeric(X%*%f$par[1:4])) }
q50<-qfit[["0.5"]]; q90<-qfit[["0.9"]]; q95<-qfit[["0.95"]]
ncruz <- sum(q90<q50)+sum(q95<q90)
cat("cruzamento:", ncruz, "\n")
if(ncruz > 0) warning("quantis condicionais nao monotonos em tau: ", ncruz,
                      " cruzamentos — verificar o ajuste antes de usar a figura")
covk<-sapply(taus,function(tt) mean(y_raw[II]<=qfit[[as.character(tt)]][II]))
fb<-R$fits[["beta 1 1"]]; mub<-plogis(as.numeric(X%*%fb$par[1:4])); shb<-exp(as.numeric(Z%*%fb$par[5:8]))
covb<-sapply(taus,function(tt) mean(y_raw[II]<=qbetamp(rep(tt,n),mub,shb)[II]))
cat("=== Calibracao Fase II (fracao abaixo do limite; alvo=tau) ===\n")
print(data.frame(tau=taus,Kuma_ajuste_tau=round(covk,3),Beta_derivado=round(covb,3)))
saveRDS(list(taus=taus,covk=covk,covb=covb),"calib_result.rds")
ff<-readRDS("farrington_result.rds"); cs<-readRDS("cusum_result.rds")
anom_cusum<-cs$flag; anom_farr<-ff$alarms_idx

## FIG 1: carta de controle por regressao quantilica
cairo_pdf("Fig_carta.pdf", width=9, height=4.2, family="DejaVu Sans")
par(mar=c(3.6,3.8,1.2,1.0), mgp=c(2.3,0.7,0), cex.axis=0.9, las=1)
plot(tnum, y_raw, type="h", col="grey72", lwd=1, xlab="Ano", ylab="Incidência de DDA",
     ylim=c(0,max(y_raw)*1.02), xaxs="i")
polygon(c(tnum,rev(tnum)), c(q95,rep(0,n)), col=rgb(0,0,0,0.05), border=NA)
abline(v=tnum[phaseI], lty=3, col="grey40")
lines(tnum, q50, col="grey45", lwd=1.2)
lines(tnum, q90, col="grey20", lwd=1.1, lty=2)
lines(tnum, q95, col="black", lwd=1.6)
exc<-II[y_raw[II]>q95[II]]; points(tnum[exc], y_raw[exc], pch=19, cex=0.65)
legend("topleft", bty="n", cex=0.82,
  legend=c("Incidência semanal","Mediana marginal (τ=0.5)","Limite de controle τ=0.90","Limite de controle τ=0.95","Excede o limite"),
  lty=c(1,1,2,1,NA), pch=c(NA,NA,NA,NA,19), lwd=c(1,1.2,1.1,1.6,NA),
  col=c("grey72","grey45","grey20","black","black"))
## Rodada 50: rotulos de fase na margem superior; dentro do grafico, "Fase II"
## colidia com a barra do pico do surto (2023-W15). Legenda: "marginal" (a curva e'
## mu_t^(tau), dadas so' as covariaveis) e ponto decimal, como no texto.
mtext("Fase I ", side=3, line=0.1, at=tnum[phaseI], adj=1, cex=0.8, col="grey40")
mtext(" Fase II", side=3, line=0.1, at=tnum[phaseI], adj=0, cex=0.8, col="grey40")
dev.off()

## FIG 2 (Fig_quantis): gerada por fig_quantis_corrigida.R, nao aqui.
## Rodada 41: este bloco produzia Fig_quantis.pdf em escala LINEAR e a partir do
## ajuste puro em tau=0.90, sobrescrevendo a figura publicada — que e' em escala
## logaritmica, como declara a legenda do artigo. Removido para que run_all.R
## nao destrua a figura correta. O gerador da Fig_quantis e' fig_quantis_corrigida.R.

## FIG 3: Farrington vs CUSUM (Fase II)
cairo_pdf("Fig_deteccao.pdf", width=9, height=2.8, family="DejaVu Sans")
par(mar=c(3.6,5.0,0.8,1.0), mgp=c(2.3,0.7,0), cex.axis=0.85, las=1)
plot(range(tnum[II]), c(0.6,2.4), type="n", xlab="Ano", ylab="", yaxt="n", xaxs="i")
axis(2, at=c(1,2), labels=c("CUSUM","Farrington"), cex.axis=0.85)
abline(h=c(1,2), col="grey88")
segments(tnum[anom_cusum],0.82,tnum[anom_cusum],1.18,lwd=2.2)
segments(tnum[anom_farr],1.82,tnum[anom_farr],2.18,lwd=2.2)
dev.off()
cat("figuras ok\n")
