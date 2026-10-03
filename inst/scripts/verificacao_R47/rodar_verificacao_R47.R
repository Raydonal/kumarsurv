## rodar_verificacao_R47.R -- executa os scripts de auditoria da Rodada 47, na ordem das
## dependencias, gravando a transcricao de cada um em resultados/verificacao_R47/.
## Rodar a partir da raiz do pacote:  LC_ALL=C.utf8 Rscript scripts/verificacao_R47/rodar_verificacao_R47.R
## Tempo: cerca de 15 min (R47_ep_varin.R, que refaz 600 hessianas, responde por quase a metade).
raiz <- normalizePath(Sys.getenv("PACOTE_VALIDACAO", "."))
setwd(file.path(raiz, "scripts"))
source("copula_ml.R")
dir_out <- file.path(DIR_RESULTADOS, "verificacao_R47"); dir.create(dir_out, showWarnings = FALSE)
ordem <- c("R47_teste_kalman.R", "R47_teste_hessiana.R", "R47_otimo_global.R", "R47_perfil_ar1_global.R",
           "R47_vero_exata.R", "R47_global_exata.R", "R47_valida_lib.R", "R47_periodo_modos.R",
           "R47_faseI.R", "R47_calib_controle.R", "R47_calib_robustez.R", "R47_ep_varin.R", "R47_tab_varin.R")
for (s in ordem) {
  t0 <- Sys.time(); cat(">>", s, "\n")
  zz <- file(file.path(dir_out, sub("\\.R$", ".txt", s)), open = "wt"); sink(zz, split = TRUE)
  r <- try(source(file.path("verificacao_R47", s), local = new.env(), echo = FALSE), silent = TRUE)
  sink(); close(zz)
  cat(if (inherits(r, "try-error")) paste("   FALHOU:", conditionMessage(attr(r, "condition"))) else "   ok",
      sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
