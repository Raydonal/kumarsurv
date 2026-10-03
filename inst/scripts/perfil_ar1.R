###############################################################################
## INTERVALO DE PERFIL DE VEROSSIMILHANCA PARA O AR(1), COM RESTRICAO DE
## INVERTIBILIDADE. O ES da hessiana nao e' reportavel (fronteira de
## estacionariedade). O perfil fixa phi, reotimiza o resto e forma o IC 95%.
## IMPORTANTE: perto de phi->1 surge um otimo ESPURIO com MA nao-invertivel
## (|MA|>1), em que AR e MA quase se cancelam; esses pontos sao excluidos por
## violarem a invertibilidade, condicao necessaria para identificacao do
## ARMA(1,1). O perfil valido considera apenas ajustes com |MA|<1.
###############################################################################
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI
I <- seq_len(phaseI); f <- R$fits[["kuma 1 1"]]
kx<-ncol(X); kz<-ncol(Z); ar_idx<-kx+kz+1; ma_idx<-ar_idx+1
phi_hat <- f$par[ar_idx]; l_max <- f$loglik
cat(sprintf("AR(1) estimado: %.4f (MA: %.4f) | loglik max: %.3f\n",
            phi_hat, f$par[ma_idx], l_max))

perfil <- function(phi0){
  nll_free <- function(pf){
    par <- f$par; par[-ar_idx] <- pf; par[ar_idx] <- phi0
    nll_copula(par, y=y[I], X=X[I,], Z=Z[I,], p=1, q=1, family="kuma", tau=0.5)
  }
  best_ll <- -Inf; best_ma <- NA_real_
  for(pert in list(f$par[-ar_idx], f$par[-ar_idx]+0.03, f$par[-ar_idx]-0.03)){
    o <- try(optim(pert, nll_free, method="BFGS", control=list(maxit=2000, reltol=1e-11)),
             silent=TRUE)
    if(inherits(o,"try-error")) next
    ma <- o$par[ar_idx]                 # MA na posicao ar_idx do vetor livre
    if(abs(ma) < 1 && -o$value > best_ll){ best_ll <- -o$value; best_ma <- ma }
  }
  c(ll=best_ll, ma=best_ma)
}

grid <- seq(0.80, 0.99, by=0.005)
ll <- numeric(length(grid)); ma <- numeric(length(grid))
for(i in seq_along(grid)){ r <- perfil(grid[i]); ll[i] <- r[1]; ma[i] <- r[2] }
crit <- qchisq(0.95, 1)/2
dev <- l_max - ll
inside <- which(is.finite(dev) & dev <= crit)
cat(sprintf("\nCriterio IC95%%: l_max - l_perfil <= %.3f\n\n", crit))
cat(sprintf("%-7s %12s %9s %7s\n","phi","loglik_perf","MA","dev"))
for(i in seq_along(grid)) if(grid[i]>=0.85 && is.finite(ll[i]))
  cat(sprintf("%-7.3f %12.3f %9.3f %7.3f%s\n", grid[i], ll[i], ma[i], dev[i],
              ifelse(dev[i]<=crit," *","")))
if(length(inside)>0){
  lo<-grid[min(inside)]; hi<-grid[max(inside)]
  cat(sprintf("\nIC 95%% de perfil (com |MA|<1) para AR(1): [%.3f, %.3f]\n", lo, hi))
}
## Rodada 41: era P=P, objeto que nunca existiu — o script abortava aqui, depois
## de imprimir o resultado correto, e nunca gravava o artefato. O perfil esta em ll.
saveRDS(list(grid=grid, lp=ll, ma=ma, l_max=l_max, phi_hat=phi_hat), "perfil_ar1.rds")
