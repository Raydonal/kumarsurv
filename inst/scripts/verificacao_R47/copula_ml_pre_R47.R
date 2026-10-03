## copula_ml.R
##
## Resolucao de caminhos (Rodada 33): todo script do pacote faz
## source("copula_ml.R"), de modo que carregar aqui o config_paths.R faz com que
## dados/, resultados/ e figuras/ sejam encontrados a partir de qualquer
## diretorio de trabalho, sem alterar os scripts.
local({
  if (exists(".PATHS_CONFIGURADOS", envir = globalenv())) return(invisible(NULL))
  for (cand in c("config_paths.R",
                 file.path("scripts", "config_paths.R"),
                 file.path("..", "scripts", "config_paths.R"))) {
    if (file.exists(cand)) { source(cand, local = globalenv()); return(invisible(NULL)) }
  }
  warning("config_paths.R nao encontrado; usando caminhos relativos ao diretorio corrente")
})

## Estimacao por maxima verossimilhanca da regressao marginal por copula
## gaussiana com dependencia ARMA, para as marginais Kumaraswamy (parametrizada
## por quantil) e beta (media-precisao). Verossimilhanca continua em forma
## fechada via a decomposicao de inovacoes (Durbin-Levinson). Sem dependencias
## de CRAN alem do R base.

## --- Kumaraswamy parametrizada por quantil (mu = tau-quantil, theta = forma) ---
.kuma_a <- function(mu, theta, tau = 0.5) log(1 - tau) / log(1 - mu^theta)
dkuma <- function(y, mu, theta, tau = 0.5, log = FALSE) {
  a  <- .kuma_a(mu, theta, tau)
  ld <- log(a) + log(theta) + (theta - 1) * log(y) + (a - 1) * log1p(-y^theta)
  if (log) ld else exp(ld)
}
pkuma <- function(y, mu, theta, tau = 0.5) {
  a <- .kuma_a(mu, theta, tau); 1 - (1 - y^theta)^a
}
qkuma <- function(u, mu, theta, tau = 0.5) {
  a <- .kuma_a(mu, theta, tau); (1 - (1 - u)^(1 / a))^(1 / theta)
}

## --- beta media-precisao (Ferrari-Cribari) ---
dbetamp <- function(y, mu, phi, log = FALSE) dbeta(y, mu * phi, (1 - mu) * phi, log = log)
pbetamp <- function(y, mu, phi) pbeta(y, mu * phi, (1 - mu) * phi)
qbetamp <- function(u, mu, phi) qbeta(u, mu * phi, (1 - mu) * phi)

## --- Durbin-Levinson (validado) ---
.dl_recursion <- function(rho) {
  n <- length(rho); phi <- vector("list", n - 1); v <- numeric(n); v[1] <- 1
  phi[[1]] <- rho[2] / v[1]; v[2] <- v[1] * (1 - phi[[1]]^2)
  if (n - 1 >= 2) for (k in 2:(n - 1)) {
    pkk <- (rho[k + 1] - sum(phi[[k - 1]] * rho[k:2])) / v[k]
    phi[[k]] <- c(phi[[k - 1]] - pkk * rev(phi[[k - 1]]), pkk)
    v[k + 1] <- v[k] * (1 - pkk^2)
  }
  list(phi = phi, v = v)
}
.one_step_moments <- function(eps, dl) {
  n <- length(eps); m <- numeric(n); s2 <- numeric(n); m[1] <- 0; s2[1] <- dl$v[1]
  for (t in 2:n) { ph <- dl$phi[[t - 1]]; m[t] <- sum(ph * eps[(t - 1):1]); s2[t] <- dl$v[t] }
  list(m = m, s2 = s2)
}

## --- log-determinante de Omega por Durbin-Levinson truncado (O(1) em n) ---
## As variancias de previsao convergem geometricamente; truncar em L e' exato
## ate precisao de maquina para L moderado.
.logdet_arma <- function(ar, ma, n, L = 80) {
  L <- min(L, n - 1)
  rho <- as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = L))
  v <- numeric(L + 1); v[1] <- 1
  phivec <- rho[2] / v[1]; v[2] <- v[1] * (1 - phivec^2)
  if (L >= 2) for (k in 2:L) {
    pkk <- (rho[k + 1] - sum(phivec * rho[k:2])) / v[k]
    phivec <- c(phivec - pkk * rev(phivec), pkk); v[k + 1] <- v[k] * (1 - pkk^2)
  }
  sum(log(v)) + (n - 1 - L) * log(v[L + 1])
}

## --- inovacoes padronizadas e forma quadratica por filtro de Kalman (O(n)) ---
.arma_kalman <- function(eps, ar, ma) {
  n <- length(eps)
  if (length(ar) == 0 && length(ma) == 0)
    return(list(quad = sum(eps^2), resid = eps))
  psi <- ARMAtoMA(ar = ar, ma = ma, lag.max = 2000); g0 <- 1 + sum(psi^2)
  mod <- makeARIMA(phi = ar, theta = ma, Delta = numeric(0))
  kr <- KalmanRun(eps * sqrt(g0), mod, update = FALSE)
  list(quad = n * kr$values[2], resid = as.numeric(kr$resid))
}

## --- log-verossimilhanca negativa da copula continua ---
nll_copula <- function(par, y, X, Z, p, q, family = "kuma", tau = 0.5) {
  kx <- ncol(X); kz <- ncol(Z)
  beta  <- par[1:kx]; gamma <- par[(kx + 1):(kx + kz)]
  ar <- if (p > 0) par[(kx + kz + 1):(kx + kz + p)] else numeric(0)
  ma <- if (q > 0) par[(kx + kz + p + 1):(kx + kz + p + q)] else numeric(0)
  mu <- plogis(as.numeric(X %*% beta))
  sh <- exp(as.numeric(Z %*% gamma))
  if (any(!is.finite(mu)) || any(mu <= 0 | mu >= 1) || any(sh <= 0)) return(1e10)
  if (family == "kuma") { fd <- dkuma(y, mu, sh, tau, log = TRUE); u <- pkuma(y, mu, sh, tau) }
  else                  { fd <- dbetamp(y, mu, sh, log = TRUE);    u <- pbetamp(y, mu, sh) }
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12); eps <- qnorm(u)
  n <- length(y)
  if (p > 0 || q > 0) {
    ak <- try(.arma_kalman(eps, ar, ma), silent = TRUE)
    ld <- try(.logdet_arma(ar, ma, n, 80), silent = TRUE)
    if (inherits(ak, "try-error") || inherits(ld, "try-error") ||
        any(!is.finite(ak$resid)) || !is.finite(ld)) return(1e10)
    quad <- ak$quad; logdet <- ld
  } else { quad <- sum(eps^2); logdet <- 0 }
  ll <- sum(fd) - 0.5 * (logdet + quad - sum(eps^2))
  if (!is.finite(ll)) return(1e10)
  -ll
}

## --- ajuste ---
fit_copula <- function(y, X, Z, p = 0, q = 0, family = "kuma", tau = 0.5) {
  kx <- ncol(X); kz <- ncol(Z)
  obj <- function(par) nll_copula(par, y = y, X = X, Z = Z, p = p, q = q,
                                  family = family, tau = tau)
  obj0 <- function(par) nll_copula(par, y = y, X = X, Z = Z, p = 0, q = 0,
                                   family = family, tau = tau)
  ## valores iniciais: modelo de independencia (BFGS; Nelder-Mead deste build e' instavel)
  st0 <- c(rep(0, kx), rep(0, kz))
  st0[1] <- if (family == "kuma") qlogis(median(y)) else qlogis(mean(y))
  ind <- optim(st0, obj0, method = "BFGS", control = list(maxit = 1000))
  start <- c(ind$par, rep(0.01, p + q))
  ## ajuste completo por BFGS, com reinicio para robustez
  fit <- optim(start, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-11))
  fit <- optim(fit$par, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  npar <- length(fit$par)
  nm <- c(paste0("mu.", colnames(X)), paste0("shape.", colnames(Z)))
  if (p > 0) nm <- c(nm, paste0("ar", 1:p)); if (q > 0) nm <- c(nm, paste0("ma", 1:q))
  names(fit$par) <- nm
  ll <- -fit$value
  se <- tryCatch({
    H <- numDeriv::hessian(obj, fit$par)
    sqrt(diag(solve(H)))
  }, error = function(e) rep(NA, npar))
  list(par = fit$par, se = se, loglik = ll, npar = npar,
       aic = -2 * ll + 2 * npar, bic = -2 * ll + log(length(y)) * npar,
       family = family, p = p, q = q, tau = tau, X = X, Z = Z, y = y,
       convergence = fit$convergence)
}

## --- residuos quantilicos preditivos de um passo ---
qresiduals <- function(fit) {
  kx <- ncol(fit$X); kz <- ncol(fit$Z); par <- fit$par
  mu <- plogis(as.numeric(fit$X %*% par[1:kx]))
  sh <- exp(as.numeric(fit$Z %*% par[(kx + 1):(kx + kz)]))
  ar <- if (fit$p > 0) par[(kx + kz + 1):(kx + kz + fit$p)] else numeric(0)
  ma <- if (fit$q > 0) par[(kx + kz + fit$p + 1):(kx + kz + fit$p + fit$q)] else numeric(0)
  u <- if (fit$family == "kuma") pkuma(fit$y, mu, sh, fit$tau) else pbetamp(fit$y, mu, sh)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12); eps <- qnorm(u)
  n <- length(fit$y)
  if (fit$p > 0 || fit$q > 0) {
    ak <- .arma_kalman(eps, ar, ma); r <- ak$resid
  } else r <- eps
  list(r = r, mu = mu, shape = sh, eps = eps)
}

## --- simulacao a partir do modelo (para validacao e bootstrap) ---
sim_copula <- function(X, Z, beta, gamma, ar = numeric(0), ma = numeric(0),
                       family = "kuma", tau = 0.5) {
  n <- nrow(X)
  if (length(ar) || length(ma)) {
    eps <- arima.sim(n = n, model = list(ar = ar, ma = ma))
    psi <- ARMAtoMA(ar = ar, ma = ma, lag.max = 1e5); eps <- eps / sqrt(1 + sum(psi^2))
  } else eps <- rnorm(n)
  mu <- plogis(as.numeric(X %*% beta)); sh <- exp(as.numeric(Z %*% gamma))
  u <- pnorm(eps)
  if (family == "kuma") qkuma(u, mu, sh, tau) else qbetamp(u, mu, sh)
}

## ---------------------------------------------------------------------------
## fit_copula_multi_tau: ajusta o modelo em varios niveis tau e verifica a
## INVARIANCIA DA LOG-VEROSSIMILHANCA MAXIMIZADA como diagnostico sistematico
## de convergencia.
##
## Fundamento: a reparametrizacao por tau (mu_t = tau-quantil condicional)
## descreve a MESMA familia de distribuicoes para qualquer tau; muda apenas a
## parametrizacao, nao o modelo. Logo, sob convergencia correta, a log-veros-
## similhanca maximizada deve ser (aproximadamente) a MESMA em todos os niveis
## tau. Um ajuste cuja log-verossimilhanca destoa dos demais nao encontrou o
## maximo global e seus quantis nao devem ser usados sem nova tentativa.
##
## Esta funcao NAO tenta corrigir automaticamente via reinicializacao em
## cascata: uma investigacao anterior mostrou que a cascata pode convergir
## para uma regiao degenerada do espaco de parametros que, apesar de log-
## verossimilhanca nominalmente mais alta, viola a propria invariancia (sinal
## inequivoco de local espurio). Em vez disso, sinaliza o problema e tenta uma
## unica reinicializacao a partir do tau de referencia mais proximo, mantendo
## sempre como resultado apenas ajustes que passam no teste de invariancia.
##
## Args:
##   taus     niveis tau a ajustar
##   tol      tolerancia absoluta em log-verossimilhanca para o teste de
##            invariancia (default: 1, escolha conservadora - ajustes que
##            diferem por mais de 1 unidade de log-verossimilhanca do valor
##            mediano sao sinalizados)
## Retorna lista nomeada por tau: $par, $loglik, $invariante (logico),
## $fonte ("padrao" ou "reinicializado"), e um atributo "diagnostico" (data
## frame) resumindo o teste para todos os niveis.
## ---------------------------------------------------------------------------
fit_copula_multi_tau <- function(y, X, Z, p = 0, q = 0, family = "kuma",
                                 taus = c(0.50, 0.75, 0.90, 0.95, 0.99),
                                 tol_ll = 1, verbose = TRUE) {
  taus <- sort(taus)
  fits <- vector("list", length(taus)); names(fits) <- as.character(taus)

  fit_um <- function(tau0, start = NULL) {
    if (is.null(start)) {
      o <- try(fit_copula(y, X, Z, p = p, q = q, family = family, tau = tau0),
               silent = TRUE)
      if (inherits(o, "try-error")) return(list(par = NA, loglik = -Inf))
      return(list(par = o$par, loglik = o$loglik))
    }
    nll <- function(par) nll_copula(par, y = y, X = X, Z = Z, p = p, q = q,
                                    family = family, tau = tau0)
    o <- try(optim(start, nll, method = "BFGS",
                   control = list(maxit = 1500, reltol = 1e-11)), silent = TRUE)
    if (inherits(o, "try-error")) return(list(par = NA, loglik = -Inf))
    list(par = o$par, loglik = -o$value)
  }

  ## 1a passada: ajuste padrao em cada nivel
  for (i in seq_along(taus)) fits[[i]] <- c(fit_um(taus[i]), fonte = "padrao")

  quantil_max <- function(par) max(plogis(as.numeric(X %*% par[1:ncol(X)])))
  qmax <- vapply(fits, function(f) if (all(is.finite(f$par))) quantil_max(f$par) else NA, numeric(1))

  ## CRITERIO ESTRUTURAL (primario): mu_t^(tau), o tau-quantil condicional, deve
  ## ser MONOTONO NAO DECRESCENTE em tau para cada t, por definicao de quantil.
  ## Isso e mais diagnostico do que a log-verossimilhanca: a superficie pode ser
  ## quase plana perto de um otimo espurio (mesma verossimilhanca, quantis
  ## implausiveis), o que a checagem de invariancia sozinha nao captura.
  mu_mat <- sapply(fits, function(f) if (all(is.finite(f$par)))
                    plogis(as.numeric(X %*% f$par[1:ncol(X)])) else rep(NA, nrow(X)))
  monotona <- function(M) all(apply(M, 1, function(r) all(diff(r) >= -1e-8), simplify = TRUE), na.rm = TRUE)
  problema <- !vapply(seq_along(taus), function(i) {
    if (i == 1 || i == length(taus)) return(TRUE)
    TRUE  # checagem completa feita abaixo, coluna a coluna
  }, logical(1))
  ## indices com violacao de monotonicidade em relacao ao vizinho anterior
  viola <- rep(FALSE, length(taus))
  for (i in 2:length(taus)) viola[i] <- any(mu_mat[, i] < mu_mat[, i - 1] - 1e-8, na.rm = TRUE)

  ## 2a passada: reinicializar os que violam monotonicidade OU cuja log-veros-
  ## similhanca destoa da mediana (ambos sao sinais de otimo espurio). Testa-se
  ## como ponto de partida CADA ajuste ainda considerado confiavel (nao apenas
  ## o tau mais proximo), pois a bacia de atracao espuria pode ser alcancada a
  ## partir de varios pontos de partida "razoaveis"; fica-se com a reinicia-
  ## lizacao de maior log-verossimilhanca dentre as que restauram a monotoni-
  ## cidade com os vizinhos ja aceitos.
  lls0 <- vapply(fits, function(f) f$loglik, numeric(1))
  ll_ref0 <- stats::median(lls0[is.finite(lls0)])
  suspeito <- viola | (abs(ll_ref0 - lls0) > tol_ll)
  ok_idx <- which(!suspeito)
  for (i in which(suspeito)) {
    if (length(ok_idx) == 0) break
    candidatos <- list()
    for (j in ok_idx) {
      novo <- fit_um(taus[i], start = fits[[j]]$par)
      if (!all(is.finite(novo$par))) next
      qnovo <- plogis(as.numeric(X %*% novo$par[1:ncol(X)]))
      mono_ok <- TRUE
      for (k in ok_idx) {
        if (taus[i] > taus[k]) mono_ok <- mono_ok && all(qnovo >= mu_mat[, k] - 1e-8)
        else                   mono_ok <- mono_ok && all(qnovo <= mu_mat[, k] + 1e-8)
      }
      if (mono_ok) candidatos[[length(candidatos) + 1]] <- list(par = novo$par,
                                                                 loglik = novo$loglik,
                                                                 q = qnovo)
    }
    if (length(candidatos) > 0) {
      melhor <- candidatos[[which.max(vapply(candidatos, `[[`, 0, "loglik"))]]
      fits[[i]] <- list(par = melhor$par, loglik = melhor$loglik, fonte = "reinicializado")
      mu_mat[, i] <- melhor$q
      ok_idx <- union(ok_idx, i)
    }
  }

  ## reavalia monotonicidade e invariancia da log-verossimilhanca apos a 2a passada
  for (i in 2:length(taus)) viola[i] <- any(mu_mat[, i] < mu_mat[, i - 1] - 1e-8, na.rm = TRUE)
  lls <- vapply(fits, function(f) f$loglik, numeric(1))
  ll_ref <- stats::median(lls[is.finite(lls)])
  ll_ok <- abs(ll_ref - lls) <= tol_ll
  aprovado <- ll_ok & !viola

  for (i in seq_along(taus)) fits[[i]]$aprovado <- aprovado[i]

  diag <- data.frame(tau = taus, loglik = round(lls, 3),
                     fonte = vapply(fits, `[[`, "", "fonte"),
                     monotono = !viola, invariante_ll = ll_ok, aprovado = aprovado)
  if (verbose) {
    cat("Diagnostico multi-tau (monotonicidade dos quantis + invariancia da log-verossimilhanca):\n")
    print(diag, row.names = FALSE)
    if (any(!aprovado))
      cat("ATENCAO: nivel(is) tau =", paste(taus[!aprovado], collapse = ", "),
          "nao aprovados mesmo apos reinicializacao; inspecionar manualmente.\n")
  }
  attr(fits, "diagnostico") <- diag
  fits
}
