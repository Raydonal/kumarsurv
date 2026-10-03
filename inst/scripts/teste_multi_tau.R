## Teste da funcao fit_copula_multi_tau nos dados da aplicacao (harmonico semestral)
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); X<-R$X; Z<-R$Z; y<-R$y; phaseI<-R$phaseI
I <- seq_len(phaseI)
fits <- fit_copula_multi_tau(y[I], X[I,], Z[I,], p=1, q=1, family="kuma",
                             taus=c(0.50,0.75,0.90,0.95,0.99), tol=1, verbose=TRUE)
saveRDS(fits, "multitau_result.rds")
