suppressMessages({library(numDeriv)})
source("copula_ml.R")
d <- read.csv("dda_platina_sem.csv", stringsAsFactors=FALSE)
n <- nrow(d)
y_raw <- d$incidence
y <- ifelse(y_raw==0, 1e-7, y_raw)            # substituicao de borda (modelo continuo)
X <- cbind(1, d$trend, d$cos26, d$sin26); Z <- X
colnames(X)<-colnames(Z)<-c("int","trend","cos26","sin26")
phaseI <- 300L
idx <- seq_len(phaseI)

berk <- function(r){ r<-r[is.finite(r)]; n<-length(r)
  nll<-function(pp){m<-pp[1];rho<-tanh(pp[2]);s2<-exp(pp[3]);z<-r-m
    -(dnorm(z[1],0,sqrt(s2/(1-rho^2)),log=TRUE)+sum(dnorm(z[-1],rho*z[-n],sqrt(s2),log=TRUE)))}
  o<-optim(c(0,0,0),nll,method="BFGS"); LR<-2*(-o$value-sum(dnorm(r,0,1,log=TRUE)))
  c(LR=LR, p=pchisq(LR,3,lower.tail=FALSE)) }

# --- ajustar Kuma e beta na Fase I, varias ordens ARMA ---
grid <- list(c(0,0),c(1,0),c(1,1),c(2,1))
res <- data.frame()
fits <- list()
for(fam in c("kuma","beta")) for(pq in grid){
  f <- try(suppressWarnings(fit_copula(y[idx],X[idx,],Z[idx,],p=pq[1],q=pq[2],family=fam)),silent=TRUE)
  if(inherits(f,"try-error")||f$convergence!=0) next
  qr <- qresiduals(f); bt <- berk(qr$r)
  res <- rbind(res, data.frame(familia=fam, p=pq[1], q=pq[2], npar=f$npar,
              loglik=round(f$loglik,2), AIC=round(f$aic,2), BIC=round(f$bic,2),
              Berk_p=round(bt["p"],3)))
  fits[[paste(fam,pq[1],pq[2])]] <- f
}
cat("=== Ajustes Fase I (n=300): Kumaraswamy vs beta ===\n")
print(res, row.names=FALSE)
saveRDS(list(res=res, fits=fits, d=d, X=X, Z=Z, y=y, phaseI=phaseI), "reais_result_sem.rds")
