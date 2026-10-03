###############################################################################
## Significância das covariáveis e especificações alternativas de tendência
## Modelo: cópula gaussiana ARMA(1,1) + marginal Kumaraswamy quantílica (Fase I)
## Rodar: Rscript analise_covariaveis.R
###############################################################################
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds")
fit <- R$fits[["kuma 1 1"]]

cat("=============================================================\n")
cat(" PARTE 1 — Significancia dos coeficientes (modelo adotado)\n")
cat(" EP por Hessiana numerica (numDeriv) da log-verossimilhanca\n")
cat("=============================================================\n")
est <- fit$par; se <- fit$se
z <- est/se; pv <- 2*pnorm(-abs(z))
nm <- if(!is.null(names(est))) names(est) else paste0("par",seq_along(est))
tab <- data.frame(parametro=nm, estimativa=round(est,4), EP=round(se,4),
                  z=round(z,2), p_valor=round(pv,4),
                  signif=ifelse(pv<0.05,"*",""), row.names=NULL)
print(tab, row.names=FALSE)
## Rodada 41: o EP do AR(1) e' NaN (fronteira de estacionariedade), o que fazia
## esta contagem sair como NA. Conta-se agora sobre os parametros avaliaveis.
cat(sprintf("\nNao significativos ao nivel de 5%%: %d de %d avaliaveis (%d sem EP)\n",
            sum(pv>=0.05, na.rm=TRUE), sum(!is.na(pv)), sum(is.na(pv))))
## Rodada 47: com o passo de hessiana corrigido o EP do AR(1) existe; o intervalo de Wald
## correspondente ultrapassa a fronteira de estacionariedade (nota da tab:estim)
iar <- grep("^ar", nm)
cat(sprintf("IC de Wald 95%% para %s: [%.3f, %.3f]\n", nm[iar], est[iar] - qnorm(0.975) * se[iar], est[iar] + qnorm(0.975) * se[iar]))

cat("\n=============================================================\n")
cat(" PARTE 2 — Especificacoes alternativas de TENDENCIA\n")
cat(" s (linear, atual) vs sqrt(s) vs log(s) vs sem tendencia\n")
cat("=============================================================\n")
d <- R$d; phaseI <- R$phaseI
y <- d$incidence[1:phaseI]
y <- ifelse(y==0, 1e-7, y)   # substituicao de borda (marginal continua)
s <- seq_len(phaseI)
trends <- list(
  "linear s"   = (s-mean(s))/100,
  "sqrt(s)"    = { v<-sqrt(s); (v-mean(v))/sd(v) },
  "log(s)"     = { v<-log(s); (v-mean(v))/sd(v) },
  "sem tend."  = NULL
)
## Rodada 50: harmonicos SEMESTRAIS (cos26, sin26), os do modelo adotado. Ate' o v13
## esta parte usava d$cos e d$sin, o harmonico ANUAL de antes da mudanca de
## periodicidade, e a comparacao de tendencias citada na Secao 6.2 referia-se a
## outro modelo. A conclusao nao muda (verificacao_R50/R50_tendencias_semestral.R
## repete a comparacao tambem com os maximos restritos a AR <= 0.97 e <= 0.996).
cs <- d$cos26[1:phaseI]; sn <- d$sin26[1:phaseI]
for(nm2 in names(trends)){
  tr <- trends[[nm2]]
  X <- if(is.null(tr)) cbind(1,cs,sn) else cbind(1,tr,cs,sn)
  f <- try(fit_copula(y, X=X, Z=matrix(1,nrow=phaseI,ncol=1),
                      family="kuma", p=1, q=1, tau=0.5), silent=TRUE)
  if(inherits(f,"try-error")){ cat(sprintf("%-10s  falhou\n", nm2)); next }
  k <- length(f$par)
  cat(sprintf("%-10s  loglik=%9.2f   npar=%2d   AIC=%10.2f\n", nm2, f$loglik, k, -2*f$loglik+2*k))
}
cat("\nNota: AIC menor = melhor. Compare para decidir a especificacao de tendencia.\n")
