###############################################################################
## 18_calib_cauda.R — SIMULACAO DEDICADA A CALIBRACAO EM CAUDA
##
## Testa a conjectura central do artigo, que ate a Rodada 32 permanecia apenas
## declarada: a de que a vantagem de calibracao do limite estimado DIRETAMENTE
## no nivel tau (marginal Kumaraswamy reparametrizada) sobre o limite DERIVADO
## de um ajuste na media (marginal beta) se acentua sob assimetria e inflacao
## de zeros -- a condicao efetivamente encontrada na aplicacao.
##
## POR QUE ESTE DESENHO E' DIFERENTE DA tab:simrec
## A tab:simrec ajusta cada marginal a dados gerados por ela mesma ("jogo em
## casa"): so' verifica se o EMV e' bem-comportado. Aqui o processo gerador
## esta' FORA DAS DUAS FAMILIAS -- e' uma logit-normal com locacao variando no
## tempo, opcionalmente inflacionada em zero -- de modo que ambas as marginais
## estao mal especificadas, como na aplicacao real. E' isso que torna a
## comparacao informativa sobre robustez.
##
## DESENHO
##   marginal geradora : logit-normal, logit(Y) ~ N(m_t, s^2)   [nao e' beta
##                       nem Kumaraswamy], com
##                       m_t = b0 + b1*t~ + b2*cos(2*pi*t/26) + b3*sin(2*pi*t/26)
##   inflacao de zeros : com probabilidade pi0, y_t = eps (= 1e-7)
##   dependencia       : copula gaussiana com ARMA(1,1) na escala latente
##   n                 : 470 (300 Fase I para ajuste + 170 Fase II para avaliar)
##   replicas          : M = 300 (padrao unificado da Secao de simulacao)
##   cenarios          : A) pi0 = 0.573  (inflacao da Fase I da aplicacao)
##                       B) pi0 = 0      (so' assimetria; controle)
##   metodos comparados: Kuma ajustada DIRETAMENTE em tau -> limite = mu_t^(tau)
##                       beta ajustada na media          -> limite = qbeta(tau,.)
##   avaliacao         : cobertura na Fase II (alvo = tau), em tau = .90/.95/.99
##
## Rodar:  LC_ALL=C.utf8 Rscript 18_calib_cauda.R
##         M_REPLICAS=5 Rscript 18_calib_cauda.R     # ensaio rapido
###############################################################################
source("copula_ml.R")

set.seed(2026)
M     <- as.integer(Sys.getenv("M_REPLICAS", "300"))
n     <- 470
nI    <- 300
taus  <- c(0.90, 0.95, 0.99)
EPS   <- 1e-7
ar1   <- 0.70    # dependencia forte, porem longe da fronteira (estabilidade)
ma1   <- -0.40

## --- locacao da marginal geradora ------------------------------------------
## Calibrada para reproduzir a escala da aplicacao: mediana da parte positiva
## da ordem de 1e-2 e cauda superior alcancando ~1e-1.
b <- c(-4.6, 0.9, -0.55, -0.15)
s_lognorm <- 1.15

tt   <- seq_len(n)
tcen <- (tt - mean(tt)) / 100
X <- cbind(1, tcen, cos(2 * pi * tt / 26), sin(2 * pi * tt / 26))
Z <- cbind(1, tcen, cos(2 * pi * tt / 26), sin(2 * pi * tt / 26))
m_t <- as.numeric(X %*% b)

I  <- seq_len(nI)
II <- (nI + 1):n

## --- gerador ----------------------------------------------------------------
## Latente ARMA(1,1) padronizado -> uniforme -> quantil da marginal geradora.
gerar <- function(pi0) {
  e <- as.numeric(arima.sim(list(ar = ar1, ma = ma1), n = n))
  e <- e / sd(e)                      # variancia unitaria (identificabilidade)
  u <- pnorm(e)
  y <- numeric(n)
  zero <- u < pi0
  y[zero] <- EPS
  if (any(!zero)) {
    v <- (u[!zero] - pi0) / (1 - pi0)
    v <- pmin(pmax(v, 1e-10), 1 - 1e-10)
    y[!zero] <- plogis(m_t[!zero] + s_lognorm * qnorm(v))
  }
  pmin(pmax(y, EPS), 1 - 1e-9)
}

## --- um cenario -------------------------------------------------------------
rodar_cenario <- function(pi0, rotulo) {
  cat(sprintf("\n=== CENARIO %s: pi0 = %.3f (%s) ===\n", rotulo, pi0,
              if (pi0 > 0) "assimetria + inflacao de zeros" else "so' assimetria"))
  cobK <- matrix(NA_real_, M, length(taus), dimnames = list(NULL, paste0("tau", taus)))
  cobB <- cobK
  convK <- convB <- 0L
  t0 <- Sys.time()

  for (i in seq_len(M)) {
    y <- gerar(pi0)

    ## --- beta ajustada na MEDIA, limite derivado por inversao ---------------
    fb <- try(fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "beta"),
              silent = TRUE)
    if (!inherits(fb, "try-error")) {
      convB <- convB + 1L
      mub  <- plogis(as.numeric(X[II, ] %*% fb$par[1:ncol(X)]))
      phib <- exp(as.numeric(Z[II, ] %*% fb$par[(ncol(X) + 1):(ncol(X) + ncol(Z))]))
      for (j in seq_along(taus))
        cobB[i, j] <- mean(y[II] <= qbetamp(taus[j], mub, phib))
    }

    ## --- Kumaraswamy ajustada DIRETAMENTE em cada tau -----------------------
    ## Partida morna a partir do nivel anterior: e' a versao economica da busca
    ## multi-partida de fit_copula_multi_tau(), necessaria porque a superficie
    ## em tau intermediario admite ponto estacionario espurio.
    ini <- NULL; okK <- TRUE
    for (j in seq_along(taus)) {
      fk <- try(if (is.null(ini))
                  fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma", tau = taus[j])
                else {
                  nll <- function(par) nll_copula(par, y = y[I], X = X[I, ], Z = Z[I, ],
                                                  p = 1, q = 1, family = "kuma", tau = taus[j])
                  o <- optim(ini, nll, method = "BFGS",
                             control = list(maxit = 1500, reltol = 1e-11))
                  list(par = o$par, loglik = -o$value)
                }, silent = TRUE)
      if (inherits(fk, "try-error")) { okK <- FALSE; break }
      ini <- fk$par
      muk <- plogis(as.numeric(X[II, ] %*% fk$par[1:ncol(X)]))
      cobK[i, j] <- mean(y[II] <= muk)
    }
    if (okK) convK <- convK + 1L

    if (i %% 25 == 0)
      cat(sprintf("  i = %d / %d  (%.1f min)\n", i, M,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  cat(sprintf("  convergidas: Kuma %d/%d, beta %d/%d | tempo %.1f min\n",
              convK, M, convB, M, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  ## --- resumo ---------------------------------------------------------------
  res <- data.frame(tau = taus,
                    cobK = colMeans(cobK, na.rm = TRUE),
                    dpK  = apply(cobK, 2, sd, na.rm = TRUE),
                    cobB = colMeans(cobB, na.rm = TRUE),
                    dpB  = apply(cobB, 2, sd, na.rm = TRUE))
  res$viesK <- res$cobK - res$tau
  res$viesB <- res$cobB - res$tau
  ## erro absoluto medio de calibracao por replica (mede dispersao, nao so' vies)
  res$eamK <- sapply(seq_along(taus), function(j) mean(abs(cobK[, j] - taus[j]), na.rm = TRUE))
  res$eamB <- sapply(seq_along(taus), function(j) mean(abs(cobB[, j] - taus[j]), na.rm = TRUE))
  ## fracao de replicas em que a Kumaraswamy fica mais proxima do alvo
  res$venceK <- sapply(seq_along(taus), function(j)
    mean(abs(cobK[, j] - taus[j]) < abs(cobB[, j] - taus[j]), na.rm = TRUE))

  cat("\n  cobertura media (alvo = tau) e vies de calibracao:\n")
  cat(sprintf("  %-6s %18s %18s %10s\n", "tau", "Kuma@tau", "beta(media)", "Kuma vence"))
  for (j in seq_along(taus))
    cat(sprintf("  %-6.2f  %.4f (%+.4f)   %.4f (%+.4f)   %6.1f%%\n",
                res$tau[j], res$cobK[j], res$viesK[j],
                res$cobB[j], res$viesB[j], 100 * res$venceK[j]))
  cat("\n  erro absoluto medio de calibracao (menor e' melhor):\n")
  for (j in seq_along(taus))
    cat(sprintf("  %-6.2f  Kuma %.4f   beta %.4f   razao beta/Kuma %.2fx\n",
                res$tau[j], res$eamK[j], res$eamB[j], res$eamB[j] / res$eamK[j]))

  list(resumo = res, cobK = cobK, cobB = cobB, pi0 = pi0,
       convK = convK, convB = convB)
}

cat("=================================================================\n")
cat(" CALIBRACAO EM CAUDA SOB MA-ESPECIFICACAO DAS DUAS MARGINAIS\n")
cat(sprintf(" gerador: logit-normal (fora das duas familias) | M = %d | n = %d\n", M, n))
cat(sprintf(" Fase I = %d, Fase II = %d | copula ARMA(1,1): ar=%.2f ma=%.2f\n",
            nI, n - nI, ar1, ma1))
cat("=================================================================\n")

## --- CENARIO C: gerador DISCRETO, fiel a' natureza dos dados ---------------
## Os cenarios A e B geram uma marginal continua, que preenche o intervalo
## entre zero e a menor incidencia positiva. A serie real NAO tem esse
## preenchimento: e' uma contagem dividida pela populacao, de modo que os
## valores possiveis sao multiplos de 10/3025 = 0.00331 e existe um VAZIO
## entre o zero e o menor positivo. Essa e' a caracteristica dominante da
## aplicacao, e nenhum gerador continuo a reproduz. O cenario C gera
## contagens binomiais negativas (sobredispersas) com media variando no tempo
## e as converte a incidencia, preservando a dependencia serial pela copula.
POP <- 3025; ESCALA <- 10
gerar_discreto <- function(size_nb, bc) {
  mu_t <- exp(as.numeric(X %*% bc))
  e <- as.numeric(arima.sim(list(ar = ar1, ma = ma1), n = n))
  e <- e / sd(e)
  u <- pmin(pmax(pnorm(e), 1e-10), 1 - 1e-10)
  cnt <- qnbinom(u, size = size_nb, mu = mu_t)
  y <- cnt * ESCALA / POP
  pmin(pmax(y, EPS), 1 - 1e-9)
}

rodar_cenario_C <- function(size_nb, bc, rotulo) {
  cat(sprintf("\n=== CENARIO %s: gerador DISCRETO (contagem NB / populacao) ===\n", rotulo))
  cobK <- matrix(NA_real_, M, length(taus), dimnames = list(NULL, paste0("tau", taus)))
  cobB <- cobK; convK <- convB <- 0L; t0 <- Sys.time()
  for (i in seq_len(M)) {
    y <- gerar_discreto(size_nb, bc)
    fb <- try(fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "beta"), silent = TRUE)
    if (!inherits(fb, "try-error")) {
      convB <- convB + 1L
      mub  <- plogis(as.numeric(X[II, ] %*% fb$par[1:ncol(X)]))
      phib <- exp(as.numeric(Z[II, ] %*% fb$par[(ncol(X) + 1):(ncol(X) + ncol(Z))]))
      for (j in seq_along(taus)) cobB[i, j] <- mean(y[II] <= qbetamp(taus[j], mub, phib))
    }
    ini <- NULL; okK <- TRUE
    for (j in seq_along(taus)) {
      fk <- try(if (is.null(ini))
                  fit_copula(y[I], X[I, ], Z[I, ], p = 1, q = 1, family = "kuma", tau = taus[j])
                else {
                  nll <- function(par) nll_copula(par, y = y[I], X = X[I, ], Z = Z[I, ],
                                                  p = 1, q = 1, family = "kuma", tau = taus[j])
                  o <- optim(ini, nll, method = "BFGS", control = list(maxit = 1500, reltol = 1e-11))
                  list(par = o$par, loglik = -o$value)
                }, silent = TRUE)
      if (inherits(fk, "try-error")) { okK <- FALSE; break }
      ini <- fk$par
      cobK[i, j] <- mean(y[II] <= plogis(as.numeric(X[II, ] %*% fk$par[1:ncol(X)])))
    }
    if (okK) convK <- convK + 1L
    if (i %% 25 == 0) cat(sprintf("  i = %d / %d  (%.1f min)\n", i, M,
                                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  cat(sprintf("  convergidas: Kuma %d/%d, beta %d/%d | tempo %.1f min\n",
              convK, M, convB, M, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res <- data.frame(tau = taus,
                    cobK = colMeans(cobK, na.rm = TRUE), cobB = colMeans(cobB, na.rm = TRUE))
  res$viesK <- res$cobK - res$tau; res$viesB <- res$cobB - res$tau
  res$eamK <- sapply(seq_along(taus), function(j) mean(abs(cobK[, j] - taus[j]), na.rm = TRUE))
  res$eamB <- sapply(seq_along(taus), function(j) mean(abs(cobB[, j] - taus[j]), na.rm = TRUE))
  res$venceK <- sapply(seq_along(taus), function(j)
    mean(abs(cobK[, j] - taus[j]) < abs(cobB[, j] - taus[j]), na.rm = TRUE))
  cat(sprintf("  %-6s %18s %18s %10s\n", "tau", "Kuma@tau", "beta(media)", "Kuma vence"))
  for (j in seq_along(taus))
    cat(sprintf("  %-6.2f  %.4f (%+.4f)   %.4f (%+.4f)   %6.1f%%\n", res$tau[j],
                res$cobK[j], res$viesK[j], res$cobB[j], res$viesB[j], 100 * res$venceK[j]))
  cat("\n  erro absoluto medio de calibracao (menor e' melhor):\n")
  for (j in seq_along(taus))
    cat(sprintf("  %-6.2f  Kuma %.4f   beta %.4f   razao beta/Kuma %.2fx\n",
                res$tau[j], res$eamK[j], res$eamB[j], res$eamB[j] / res$eamK[j]))
  list(resumo = res, cobK = cobK, cobB = cobB, convK = convK, convB = convB)
}

A <- rodar_cenario(0.573, "A")
B <- rodar_cenario(0.000, "B")
## Parametros escolhidos por ajuste a's estatisticas descritivas da serie real
## (Platina): zeros 48.9%, menor positivo 0.00331 (= 1 caso), mediana positiva
## 0.00661 (= 2 casos), maximo 0.2380, razao maximo/mediana-positiva 36x.
## Com size=0.5 e mu0=1.5 obtem-se 49.9% / 0.00331 / 0.00661 / 0.198 / 30x.
## O quantil 0.99 fica abaixo do real (0.066 contra 0.137), de modo que o
## cenario e' se algo MENOS extremo do que os dados, e portanto conservador.
bc <- c(log(1.5), 0.3, -0.55, -0.15)
C <- rodar_cenario_C(0.50, bc, "C")

saveRDS(list(A = A, B = B, C = C, M = M, n = n, nI = nI, taus = taus,
             ar1 = ar1, ma1 = ma1, b = b, s = s_lognorm, bc = bc),
        "calib_cauda_result.rds")

cat("\n=== RESUMO FINAL ===\n")
for (cen in list(list(r = A, nm = "A (inflacionado, pi0=0.573)"),
                 list(r = B, nm = "B (so' assimetria)"))) {
  cat(sprintf("\n%s\n", cen$nm))
  r <- cen$r$resumo
  for (j in seq_len(nrow(r)))
    cat(sprintf("  tau=%.2f  Kuma %+.4f  beta %+.4f  (EAM %.4f vs %.4f)\n",
                r$tau[j], r$viesK[j], r$viesB[j], r$eamK[j], r$eamB[j]))
}
cat("\nFIM\n")
