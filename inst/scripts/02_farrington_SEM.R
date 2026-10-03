suppressMessages(library(surveillance))
source("copula_ml.R")
R <- readRDS("reais_result_sem.rds")
d <- R$d; phaseI <- R$phaseI; n <- nrow(d)

# --- construir objeto sts com as contagens semanais ---
start_year <- as.integer(substr(d$epiweek[1],1,4))
start_week <- as.integer(substr(d$epiweek[1],7,8))
counts <- d$cases
sts_obj <- sts(observed = counts, start = c(start_year, start_week), frequency = 52)

# --- Farrington Flexible (Noufaily 2013) na Fase II (semanas 301:470) ---
con <- list(range = (phaseI+1):n, b = 4, w = 3, weightsThreshold = 2.58,
            pastWeeksNotIncluded = 26, thresholdMethod = "nbPlugin",
            alpha = 0.005, trend = TRUE, pThresholdTrend = 0.05,
            noPeriods = 10, limit54 = c(0,4))
ff <- try(farringtonFlexible(sts_obj, con), silent = TRUE)
if(inherits(ff,"try-error")){ cat("erro farrington:", conditionMessage(attr(ff,"condition")),"\n") } else {
  alarms_ff <- which(as.logical(alarms(ff)))            # indices dentro do range
  wk_ff <- d$epiweek[(phaseI+1):n][alarms_ff]
  cat("=== Farrington Flexible (Fase II) ===\n")
  cat("numero de alarmes:", length(alarms_ff), "\n")
  cat("semanas sinalizadas:", paste(wk_ff, collapse=", "), "\n")
  saveRDS(list(alarms_idx=alarms_ff+phaseI, weeks=wk_ff, upperbound=upperbound(ff)), "farrington_result.rds")
}

# --- CUSUM da copula (Kuma ARMA(1,1)) com residuos prospectivos, para comparar ---
X<-R$X; Z<-R$Z; y<-R$y
fitI <- R$fits[["kuma 1 1"]]
# congelar parametros da Fase I, avaliar residuos preditivos na serie completa
par<-fitI$par; kx<-ncol(X); kz<-ncol(Z)
mu<-plogis(as.numeric(X%*%par[1:kx])); sh<-exp(as.numeric(Z%*%par[(kx+1):(kx+kz)]))
ar<-par[(kx+kz+1):(kx+kz+1)]; ma<-par[(kx+kz+2)]
u<-pkuma(y,mu,sh,0.5); u<-pmin(pmax(u,1e-12),1-1e-12); eps<-qnorm(u)
# Durbin-Levinson puro-R (mesmo caminho de 08_vantagem_SEM.R), em vez de
# .arma_kalman/stats::KalmanRun: com o AR(1) proximo da fronteira de
# estacionariedade (0.948), KalmanRun mostrou-se sensivel a versao do R;
# a recursao DL e' matematicamente equivalente e nao usa codigo interno do
# pacote stats, dando resultado estavel entre ambientes (ver RESPOSTAS_COAUTOR.md,
# Rodada 30).
rho <- as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = length(eps) - 1))
mom <- .one_step_moments(eps, .dl_recursion(rho))
r_full <- (eps - mom$m) / sqrt(mom$s2)
# CUSUM unilateral com reset na Fase II
k<-0.5; cp<-numeric(n)
for(i in (phaseI+1):n) cp[i]<-max(0, r_full[i]-k+cp[i-1])
flag_cusum <- which(cp[(phaseI+1):n] > 4) + phaseI
cat("\n=== CUSUM copula Kuma (h=4, Fase II) ===\n")
cat("numero de alarmes:", length(flag_cusum), "\n")
cat("semanas sinalizadas:", paste(d$epiweek[flag_cusum], collapse=", "), "\n")
saveRDS(list(r_full=r_full, cp=cp, flag=flag_cusum), "cusum_result.rds")
