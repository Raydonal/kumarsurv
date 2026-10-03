## O efeito real da estrutura de correlacao: a banda CONDICIONAL.
## Ignorar a dependencia produz banda com largura errada (usa sd marginal=1
## em vez do sd de inovacao de 1 passo), o que custa PODER de deteccao.
source("copula_ml.R"); set.seed(21)
n<-400; R<-400
estr<-list("independente"=list(ar=numeric(0),ma=numeric(0)),
           "AR(1) 0,3"=list(ar=0.3,ma=numeric(0)),
           "AR(1) 0,7"=list(ar=0.7,ma=numeric(0)),
           "ARMA(1,1) 0,6/0,4"=list(ar=0.6,ma=0.4))
res<-data.frame()
for(nm in names(estr)){ g<-estr[[nm]]
  # sd da inovacao de 1 passo (assintotico)
  v1<-if(length(g$ar)||length(g$ma)){rho<-as.numeric(ARMAacf(ar=g$ar,ma=g$ma,lag.max=60));.dl_recursion(rho)$v[60]} else 1
  sd1<-sqrt(v1)
  # largura relativa da banda de 95%: ignorando (sd=1) vs modelando (sd=sd1)
  larg_ig<-2*qnorm(0.975)*1; larg_ok<-2*qnorm(0.975)*sd1
  # poder: prob. de detectar um deslocamento de +1.5 na escala latente
  pot_ig<-1-pnorm(qnorm(0.95)*1,      mean=1.5, sd=sd1)   # limite calibrado errado (sd=1), erro real sd1
  pot_ok<-1-pnorm(qnorm(0.95)*sd1,    mean=1.5, sd=sd1)   # limite correto
  # taxa de alarme falso do limite calibrado errado
  fa_ig<-1-pnorm(qnorm(0.95)*1, mean=0, sd=sd1)
  res<-rbind(res,data.frame(estrutura=nm, sd_inovacao=round(sd1,3),
    largura_ignora=round(larg_ig,2), largura_correta=round(larg_ok,2),
    alarme_falso_ignora=round(fa_ig,4), alarme_falso_correto=0.05,
    poder_ignora=round(pot_ig,3), poder_correto=round(pot_ok,3)))
}
cat("=== Efeito da estrutura de correlacao da copula sobre a banda e a carta ===\n")
cat("(banda de 95%; deslocamento de 1,5 na escala latente)\n\n")
print(res,row.names=FALSE)
saveRDS(res,"corr_result.rds")
