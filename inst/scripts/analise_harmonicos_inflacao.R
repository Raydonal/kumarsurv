###############################################################################
## (A) Harmônicos: período ANUAL vs períodos INFERIORES a um ano
##     (decide a questão "sazonalidade vs ciclo")
## (B) Modelo INFLACIONADO em zero COM dependência ARMA via cópula
##     (hurdle: P(Y=0)=pi ; Y|Y>0 ~ marginal + cópula gaussiana ARMA)
## Rodar: Rscript analise_harmonicos_inflacao.R
###############################################################################
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); d <- R$d; phaseI <- R$phaseI
y_raw <- d$incidence[1:phaseI]
y <- ifelse(y_raw==0, 1e-7, y_raw)
s <- seq_len(phaseI)
tr <- (s-mean(s))/100

cat("=============================================================\n")
cat(" (A) PERIODICIDADE DOS HARMONICOS\n")
cat(" H0 pratico: qual periodo ajusta melhor? 52 sem (anual, = CICLO)\n")
cat(" vs 26 (semestral), 13 (trimestral), 4.33 (mensal) = SAZONAIS\n")
cat("=============================================================\n")
per <- c("52 (anual)"=52, "26 (semestral)"=26, "13 (trimestral)"=13, "4.33 (mensal)"=365.25/12/7)
res <- data.frame()
for(nm in names(per)){
  P <- per[[nm]]
  X <- cbind(1, tr, cos(2*pi*s/P), sin(2*pi*s/P))
  f <- try(fit_copula(y, X=X, Z=matrix(1,phaseI,1), family="kuma", p=1, q=1, tau=0.5), silent=TRUE)
  if(inherits(f,"try-error")){ cat(sprintf("%-15s falhou\n",nm)); next }
  k <- length(f$par)
  res <- rbind(res, data.frame(periodo=nm, loglik=round(f$loglik,2), npar=k,
                               AIC=round(-2*f$loglik+2*k,2)))
}
# sem harmonicos, para referencia
f0 <- try(fit_copula(y, X=cbind(1,tr), Z=matrix(1,phaseI,1), family="kuma", p=1, q=1, tau=0.5), silent=TRUE)
if(!inherits(f0,"try-error")){
  k0<-length(f0$par)
  res <- rbind(res, data.frame(periodo="sem harmonicos", loglik=round(f0$loglik,2),
                               npar=k0, AIC=round(-2*f0$loglik+2*k0,2)))
}
print(res, row.names=FALSE)
cat("\nLeitura: se nenhum periodo <1 ano melhora o ajuste e o modelo sem\n")
cat("harmonicos nao e' pior, nao ha' evidencia de SAZONALIDADE nesta serie.\n")

cat("\n=============================================================\n")
cat(" (B) MODELO INFLACIONADO **COM DEPENDENCIA ARMA** (hurdle+copula)\n")
cat(" ll_total = [n0*log(pi) + npos*log(1-pi)] + ll_copula(positivos)\n")
cat("=============================================================\n")
is0 <- y_raw==0; ypos <- y_raw[!is0]; npos <- length(ypos); n0 <- sum(is0); n <- phaseI
pihat <- n0/n
ll_infl <- n0*log(pihat) + npos*log(1-pihat)
spos <- s[!is0]; trpos <- (spos-mean(spos))/100
Xp <- cbind(1, trpos, cos(2*pi*spos/52), sin(2*pi*spos/52))
cat(sprintf("pi_hat=%.4f | n0=%d | npos=%d\n", pihat, n0, npos))
cat("Obs.: os positivos formam subsequencia irregularmente espacada; a estrutura\n")
cat("ARMA e' ajustada sobre essa subsequencia (aproximacao declarada).\n\n")
out <- data.frame()
for(fam in c("kuma","beta")){
  f <- try(fit_copula(ypos, X=Xp, Z=matrix(1,npos,1), family=fam, p=1, q=1, tau=0.5), silent=TRUE)
  if(inherits(f,"try-error")){ cat(fam,"falhou\n"); next }
  k <- length(f$par) + 1                    # +1 pelo pi
  ll <- ll_infl + f$loglik
  out <- rbind(out, data.frame(marginal=fam, ll_copula_pos=round(f$loglik,2),
               ll_total=round(ll,2), npar=k, AIC=round(-2*ll+2*k,2)))
}
print(out, row.names=FALSE)
if(nrow(out)==2) cat(sprintf("\nDiferenca AIC (beta - kuma) = %.2f (negativo => beta melhor)\n",
                             out$AIC[out$marginal=="beta"]-out$AIC[out$marginal=="kuma"]))
