## R50_tau090_pontos.R -- Rodada 50 (verificacao; rodar por rodar_verificacao_R50.R).
## Os tres pontos a que chegam os ajustes em tau = 0.90 (Kumaraswamy, Fase I, modelo da
## Tabela 2 do artigo), com log-verossimilhanca, componente autorregressivo e faixa do
## quantil de nivel 0.90:
##   (a) fit_copula_multi_tau(), o reportado na Tabela 3 (reinicio a partir de tau = 0.75);
##   (b) reinicio a partir do ajuste em tau = 0.95, o segundo ponto estacionario da nota;
##   (c) fit_copula() com partida padrao, que vai 'a bacia da fronteira de estacionariedade.
## Sustenta a nota da Tabela 3 e o paragrafo "O nivel tau = 0.90" da Secao S6.6.
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds"); X <- R$X; Z <- R$Z; y <- R$y; I <- seq_len(R$phaseI)
ymax <- max(R$d$incidence)
fm <- fit_copula_multi_tau(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma",
                           taus = c(0.50, 0.75, 0.90, 0.95, 0.99), verbose = FALSE)
a <- fm[["0.9"]]
p95 <- fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma", tau = 0.95)
nll <- function(par) nll_copula(par, y = y[I], X = X[I, ], Z = Z[I, ], p = 1, q = 1, family = "kuma", tau = 0.90)
o <- optim(p95$par, nll, method = "BFGS", control = list(maxit = 1500, reltol = 1e-11))
b <- list(par = o$par, loglik = -o$value)
c0 <- fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma", tau = 0.90)
lin <- function(nm, f) {
  q <- plogis(as.numeric(X %*% f$par[1:ncol(X)]))
  cat(sprintf("%-44s loglik %9.3f | AR %.4f | intercepto %7.3f | quantil 0.90 em [%.4f, %.4f]%s\n",
              nm, f$loglik, f$par[9], f$par[1], min(q), max(q),
              if (min(q) > ymax) "  (acima do maximo observado em todas as semanas)" else ""))
}
cat(sprintf("=== tau = 0.90: pontos estacionarios (maximo observado da serie: %.4f) ===\n", ymax))
lin("(a) multi-tau, reportado na Tabela 3", a)
lin("(b) reinicio a partir de tau = 0.95", b)
lin("(c) fit_copula, partida padrao", c0)
