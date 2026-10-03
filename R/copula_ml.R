## Maximum-likelihood estimation of the Gaussian-copula marginal regression with
## Kumaraswamy (quantile-parametrized) and beta (mean-precision) margins. The
## likelihood stays closed-form via the innovations (Durbin-Levinson) decomposition.
## No dependencies beyond base R (numDeriv is used only for standard errors, if
## installed). This file is copula_ml.R from the paper's validation package
## (inst/scripts/), stripped of the runtime path-resolution shim that script
## sourcing needs (config_paths.R) but a package does not, and documented with
## roxygen2 for export.

.kuma_a <- function(mu, theta, tau = 0.5) log(1 - tau) / log(1 - mu^theta)

#' Kumaraswamy density, quantile-parametrized
#'
#' Density of the Kumaraswamy distribution reparametrized so that `mu` is its
#' `tau`-quantile, `mu = F^{-1}(tau)`, rather than one of the two standard shape
#' parameters. For `tau = 0.5`, `mu` is the median, the parametrization of
#' Mitnik and Baek (2013).
#'
#' @param y Numeric vector in `(0, 1)`.
#' @param mu Numeric vector, the `tau`-quantile, in `(0, 1)`.
#' @param theta Numeric vector, the shape parameter, `> 0`.
#' @param tau Quantile level, in `(0, 1)`. Default `0.5` (median).
#' @param log Logical; if `TRUE`, the log-density is returned.
#' @return Numeric vector of (log-)densities.
#' @seealso [pkuma()], [qkuma()]
#' @import stats
#' @export
dkuma <- function(y, mu, theta, tau = 0.5, log = FALSE) {
  a  <- .kuma_a(mu, theta, tau)
  ld <- log(a) + log(theta) + (theta - 1) * log(y) + (a - 1) * log1p(-y^theta)
  if (log) ld else exp(ld)
}

#' Kumaraswamy distribution function, quantile-parametrized
#'
#' @inheritParams dkuma
#' @return Numeric vector, `F(y)`.
#' @seealso [dkuma()], [qkuma()]
#' @export
pkuma <- function(y, mu, theta, tau = 0.5) {
  a <- .kuma_a(mu, theta, tau); 1 - (1 - y^theta)^a
}

#' Kumaraswamy quantile function, quantile-parametrized
#'
#' Closed form: no numerical inversion is needed, unlike the beta distribution.
#' This is the property the paper's copula construction, predictive quantile
#' residuals and prediction bands rely on throughout.
#'
#' @param u Numeric vector of probabilities in `(0, 1)`.
#' @inheritParams dkuma
#' @return Numeric vector, `F^{-1}(u)`.
#' @seealso [dkuma()], [pkuma()]
#' @export
qkuma <- function(u, mu, theta, tau = 0.5) {
  a <- .kuma_a(mu, theta, tau); (1 - (1 - u)^(1 / a))^(1 / theta)
}

#' Beta density, mean-precision parametrization
#'
#' The mean-precision beta regression parametrization of Ferrari and
#' Cribari-Neto (2004), used throughout as the comparison marginal distribution.
#'
#' @param y Numeric vector in `(0, 1)`.
#' @param mu Numeric vector, the mean, in `(0, 1)`.
#' @param phi Numeric vector, the precision parameter, `> 0`.
#' @param log Logical; if `TRUE`, the log-density is returned.
#' @return Numeric vector of (log-)densities.
#' @seealso [pbetamp()], [qbetamp()]
#' @export
dbetamp <- function(y, mu, phi, log = FALSE) dbeta(y, mu * phi, (1 - mu) * phi, log = log)

#' Beta distribution function, mean-precision parametrization
#' @inheritParams dbetamp
#' @return Numeric vector, `F(y)`.
#' @export
pbetamp <- function(y, mu, phi) pbeta(y, mu * phi, (1 - mu) * phi)

#' Beta quantile function, mean-precision parametrization
#'
#' Unlike [qkuma()], this requires the numerical inversion built into
#' [stats::qbeta()]: the beta distribution function has no closed form.
#'
#' @param u Numeric vector of probabilities in `(0, 1)`.
#' @inheritParams dbetamp
#' @return Numeric vector, `F^{-1}(u)`.
#' @export
qbetamp <- function(u, mu, phi) qbeta(u, mu * phi, (1 - mu) * phi)

## --- Durbin-Levinson recursion (innovations algorithm) ---
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

## --- log-determinant of the ARMA correlation matrix, by Durbin-Levinson ---
## Runs until the prediction variance settles (|v_k - v_{k-1}| <= 1e-14 v_k),
## not to a fixed truncation length: near the stationarity boundary a fixed L
## underestimates the log-determinant (see the paper's pitfall 16).
.logdet_arma <- function(ar, ma, n, L = 80) {
  Lmax <- n - 1
  rho <- as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = Lmax))
  v <- numeric(Lmax + 1); v[1] <- 1
  phivec <- rho[2] / v[1]; v[2] <- v[1] * (1 - phivec^2); k <- 1
  while (k < Lmax) {
    k <- k + 1
    pkk <- (rho[k + 1] - sum(phivec * rho[k:2])) / v[k]
    phivec <- c(phivec - pkk * rev(phivec), pkk); v[k + 1] <- v[k] * (1 - pkk^2)
    if (k >= L && abs(v[k + 1] - v[k]) <= 1e-14 * v[k + 1]) break
  }
  sum(log(v[1:(k + 1)])) + (n - 1 - k) * log(v[k + 1])
}

## --- exact variance of an ARMA process with unit-variance innovations ---
## gamma(0) = 1 + sum(psi_j^2), solved exactly via the linear system behind
## stats::ARMAacf before normalization: matches a direct sum of psi-weights to
## 1e-12 relative error even at AR = 0.9999, at O(p^3) instead of needing
## millions of summed terms near the stationarity boundary.
.arma_var <- function(ar, ma) {
  p <- length(ar); q <- length(ma)
  if (p == 0) return(1 + sum(ma^2))
  if (any(Mod(polyroot(c(1, -ar))) <= 1)) stop("non-stationary autoregressive part")
  r <- max(p, q + 1)
  if (r == 1) return(1 / (1 - ar^2))
  if (r > p) { ar <- c(ar, rep(0, r - p)); p <- r }
  p1 <- p + 1L; p2.1 <- p + p1
  A <- matrix(0, p1, p2.1)
  ind <- seq_len(p1)
  ind <- as.matrix(expand.grid(ind, ind))[, 2L:1L]
  ind[, 2] <- ind[, 1L] + ind[, 2L] - 1L
  A[ind] <- c(1, -ar)
  A[, 1L:p] <- A[, 1L:p] + A[, p2.1:(p + 2L)]
  rhs <- c(1, rep(0, p))
  if (q > 0) {
    psi <- c(1, ARMAtoMA(ar, ma, q))
    theta <- c(1, ma, rep(0, q + 1L))
    for (k in 1L + 0:q) rhs[k] <- sum(psi * theta[k + 0:q])
  }
  ind <- p1:1
  solve(A[ind, ind], rhs)[1L]
}

## --- standardized innovations and quadratic form via the Kalman filter, O(n) ---
.arma_kalman <- function(eps, ar, ma) {
  n <- length(eps)
  if (length(ar) == 0 && length(ma) == 0)
    return(list(quad = sum(eps^2), resid = eps))
  g0 <- .arma_var(ar, ma)
  mod <- makeARIMA(phi = ar, theta = ma, Delta = numeric(0))
  kr <- KalmanRun(eps * sqrt(g0), mod, update = FALSE)
  list(quad = n * kr$values[2], resid = as.numeric(kr$resid))
}

#' Negative log-likelihood of the Gaussian-copula marginal regression
#'
#' The likelihood underlying [fit_copula()]: continuous Kumaraswamy or beta
#' margins, serial dependence induced by a Gaussian copula with ARMA(`p`,`q`)
#' structure on the latent (probit-transformed) scale. Evaluated in closed form
#' via the innovations decomposition (no explicit construction or inversion of
#' the `n x n` correlation matrix).
#'
#' @param par Numeric vector of parameters, in the order: location-submodel
#'   coefficients (length `ncol(X)`), shape-submodel coefficients (length
#'   `ncol(Z)`), AR coefficients (length `p`), MA coefficients (length `q`).
#' @param y Numeric response vector in `(0, 1)`.
#' @param X Location-submodel design matrix.
#' @param Z Shape-submodel design matrix.
#' @param p,q Non-negative integers, the ARMA orders.
#' @param family `"kuma"` (default) or `"beta"`.
#' @param tau Quantile level for the Kumaraswamy margin (ignored if `family = "beta"`).
#' @return The negative log-likelihood, a single number (`1e10` as a penalty for
#'   parameter values outside the admissible region, e.g. a non-stationary
#'   ARMA polynomial).
#' @seealso [fit_copula()]
#' @export
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

#' Fit the Gaussian-copula marginal regression by maximum likelihood
#'
#' Fits the model of [nll_copula()] by BFGS (Nelder-Mead is unstable on this
#' likelihood surface; see the paper's pitfall 6), starting from the
#' independence-model fit, then computes standard errors from the numerical
#' Hessian with a relative step of `1e-3` (the `numDeriv` default, `d = 0.1`,
#' can push an autoregressive coefficient outside the stationary region near
#' the boundary; see the paper's pitfall 16).
#'
#' **At `tau` other than the median, prefer [fit_copula_multi_tau()]**: the
#' likelihood surface can have more than one interior stationary point at a
#' given `tau`, with different values, and a bare single-start `fit_copula()`
#' call is not guaranteed to land on the one with the higher likelihood (see
#' the paper's pitfall 1).
#'
#' @inheritParams nll_copula
#' @return A list with components `par` (named numeric vector), `se` (standard
#'   errors, `NA` if the Hessian is not invertible), `loglik`, `npar`, `aic`,
#'   `bic`, `family`, `p`, `q`, `tau`, the inputs `X`, `Z`, `y`, and `convergence`
#'   (as returned by [stats::optim()]).
#' @seealso [fit_copula_multi_tau()], [qresiduals()], [sim_copula()]
#' @export
fit_copula <- function(y, X, Z, p = 0, q = 0, family = "kuma", tau = 0.5) {
  kx <- ncol(X); kz <- ncol(Z)
  obj <- function(par) nll_copula(par, y = y, X = X, Z = Z, p = p, q = q,
                                  family = family, tau = tau)
  obj0 <- function(par) nll_copula(par, y = y, X = X, Z = Z, p = 0, q = 0,
                                   family = family, tau = tau)
  st0 <- c(rep(0, kx), rep(0, kz))
  st0[1] <- if (family == "kuma") qlogis(median(y)) else qlogis(mean(y))
  ind <- optim(st0, obj0, method = "BFGS", control = list(maxit = 1000))
  start <- c(ind$par, rep(0.01, p + q))
  fit <- optim(start, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-11))
  fit <- optim(fit$par, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  npar <- length(fit$par)
  nm <- c(paste0("mu.", colnames(X)), paste0("shape.", colnames(Z)))
  if (p > 0) nm <- c(nm, paste0("ar", 1:p)); if (q > 0) nm <- c(nm, paste0("ma", 1:q))
  names(fit$par) <- nm
  ll <- -fit$value
  se <- tryCatch({
    H <- numDeriv::hessian(obj, fit$par, method.args = list(d = 1e-3))
    sqrt(diag(solve(H)))
  }, error = function(e) rep(NA, npar))
  list(par = fit$par, se = se, loglik = ll, npar = npar,
       aic = -2 * ll + 2 * npar, bic = -2 * ll + log(length(y)) * npar,
       family = family, p = p, q = q, tau = tau, X = X, Z = Z, y = y,
       convergence = fit$convergence)
}

#' Predictive quantile residuals
#'
#' One-step-ahead predictive quantile residuals (Masarotto and Varin, 2012,
#' 2017) from a model fitted by [fit_copula()]: the standardized innovations
#' evaluated at the maximum-likelihood estimates. Under correct specification
#' and known parameters, these are i.i.d. standard normal; see the paper's
#' Proposition 1 for the finite-sample statement.
#'
#' @param fit A list as returned by [fit_copula()].
#' @return A list with components `r` (the residuals), `mu`, `shape` (the
#'   fitted marginal parameters) and `eps` (the untransformed latent scale).
#' @seealso [fit_copula()]
#' @export
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

#' Simulate from the Gaussian-copula marginal regression
#'
#' Used for Monte Carlo validation and for the parametric bootstrap behind the
#' paper's goodness-of-fit tests and prediction bands.
#'
#' @param X Location-submodel design matrix.
#' @param Z Shape-submodel design matrix.
#' @param beta,gamma Numeric vectors of location/shape-submodel coefficients.
#' @param ar,ma Numeric vectors of ARMA coefficients (empty for independence).
#' @param family `"kuma"` (default) or `"beta"`.
#' @param tau Quantile level for the Kumaraswamy margin (ignored if `family = "beta"`).
#' @return Numeric vector of length `nrow(X)`, a simulated series in `(0, 1)`.
#' @export
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

#' Fit the model at several quantile levels, with a monotonicity/invariance check
#'
#' Fits [fit_copula()] at each level in `taus` and runs two systematic
#' diagnostics the paper relies on to avoid a known failure mode (pitfall 1:
#' a single-start fit at `tau = 0.90` can land on a spurious stationary point
#' whose fitted quantile sits above every observation):
#'
#' 1. **Log-likelihood invariance.** The `tau`-reparametrization describes the
#'    *same* distribution family for every `tau`; only the parametrization
#'    changes, not the model. Under correct convergence, the maximized
#'    log-likelihood should therefore be (approximately) the same across all
#'    levels. A fit whose log-likelihood departs from the others by more than
#'    `tol_ll` has not found the global maximum.
#' 2. **Quantile monotonicity.** `mu_t^{(tau)}`, the conditional `tau`-quantile,
#'    must be non-decreasing in `tau` for every `t`, by definition of a
#'    quantile. This catches spurious optima the invariance check alone can
#'    miss (the likelihood surface can be nearly flat near a spurious optimum).
#'
#' Fits flagged by either check are retried once from the parameter estimates
#' of each already-accepted level (not only the nearest one, since the
#' spurious basin can be reached from several "reasonable" starting points),
#' keeping the restart with the highest log-likelihood among those that
#' restore monotonicity with the levels already accepted. This function does
#' **not** cascade restarts further: an earlier investigation found that a
#' restart cascade can converge to a degenerate region of the parameter space
#' that, despite a nominally higher log-likelihood, violates the invariance
#' property itself (an unambiguous sign of a spurious local optimum) --
#' instead, unresolved cases are flagged for manual inspection.
#'
#' @inheritParams fit_copula
#' @param taus Numeric vector of quantile levels to fit.
#' @param tol_ll Absolute log-likelihood tolerance for the invariance check
#'   (default `1`; fits whose log-likelihood differs from the median by more
#'   than this are flagged).
#' @param verbose Logical; print the diagnostic table.
#' @return A list named by `tau`, each element with components `par`,
#'   `loglik`, `fonte` (`"padrao"` or `"reinicializado"`: standard fit or
#'   restarted) and `aprovado` (logical, passed both checks), with a
#'   `"diagnostico"` attribute (a data frame) summarizing the checks for every
#'   level.
#' @seealso [fit_copula()]
#' @export
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

  for (i in seq_along(taus)) fits[[i]] <- c(fit_um(taus[i]), fonte = "padrao")

  quantil_max <- function(par) max(plogis(as.numeric(X %*% par[1:ncol(X)])))
  qmax <- vapply(fits, function(f) if (all(is.finite(f$par))) quantil_max(f$par) else NA, numeric(1))

  mu_mat <- sapply(fits, function(f) if (all(is.finite(f$par)))
                    plogis(as.numeric(X %*% f$par[1:ncol(X)])) else rep(NA, nrow(X)))
  viola <- rep(FALSE, length(taus))
  for (i in 2:length(taus)) viola[i] <- any(mu_mat[, i] < mu_mat[, i - 1] - 1e-8, na.rm = TRUE)

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
    cat("Multi-tau diagnostic (quantile monotonicity + log-likelihood invariance):\n")
    print(diag, row.names = FALSE)
    if (any(!aprovado))
      cat("WARNING: level(s) tau =", paste(taus[!aprovado], collapse = ", "),
          "not approved even after restart; inspect manually.\n")
  }
  attr(fits, "diagnostico") <- diag
  fits
}
