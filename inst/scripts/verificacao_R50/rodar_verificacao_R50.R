## rodar_verificacao_R50.R -- executa as duas verificacoes da Rodada 50, gravando a
## transcricao de cada uma em resultados/verificacao_R50/.
## Rodar a partir da raiz do pacote:  LC_ALL=C.utf8 Rscript scripts/verificacao_R50/rodar_verificacao_R50.R
## Tempo: cerca de 7 min (R50_tendencias_semestral.R responde por quase todo).
raiz <- normalizePath(Sys.getenv("PACOTE_VALIDACAO", "."))
setwd(file.path(raiz, "scripts"))
source("copula_ml.R")
dir_out <- file.path(DIR_RESULTADOS, "verificacao_R50"); dir.create(dir_out, showWarnings = FALSE)
for (s in c("R50_tau090_pontos.R", "R50_tendencias_semestral.R")) {
  t0 <- Sys.time(); cat(">>", s, "\n")
  zz <- file(file.path(dir_out, sub("\\.R$", ".txt", s)), open = "wt"); sink(zz, split = TRUE)
  r <- try(source(file.path("verificacao_R50", s), local = new.env(), echo = FALSE), silent = TRUE)
  sink(); close(zz)
  cat(if (inherits(r, "try-error")) paste("   FALHOU:", conditionMessage(attr(r, "condition"))) else "   ok",
      sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
