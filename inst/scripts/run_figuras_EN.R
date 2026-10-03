## Regenerates, with English labels, the 8 figures actually included in
## v6_rodada50/main_v14.tex and supplementary_v14.tex, writing PDFs straight
## into v6_rodada50/ImgsIng/. Run from PACOTE_VALIDACAO/scripts/:
##   LC_ALL=C.utf8 LANG=C.utf8 Rscript run_figuras_EN.R
source("copula_ml.R")  # loads config_paths.R, defines DIR_SCRIPTS etc.
setwd(DIR_SCRIPTS)

DIR_FIGURAS_EN <<- normalizePath(file.path(RAIZ, "..", "v6_rodada50", "ImgsIng"), mustWork = FALSE)
if (!dir.exists(DIR_FIGURAS_EN)) dir.create(DIR_FIGURAS_EN, recursive = TRUE)
cat("[run_figuras_EN] output dir:", DIR_FIGURAS_EN, "\n")

scripts_en <- c(
  "fig_densidades_EN.R",       # Fig_densidades (supplementary)
  "fig3_horizontal_EN.R",      # Fig3 (main)
  "fig4_cusum_SEM_EN.R",       # Fig4 (main)
  "06_figuras_SEM_EN.R",       # Fig_carta (main)
  "fig_quantis_corrigida_EN.R",# Fig_quantis (supplementary)
  "07_fig_deteccao_SEM_EN.R",  # Fig_deteccao (main)
  "11_arl_EN.R",                # Fig_arl (supplementary)
  "10_banda_simulada_EN.R"      # Fig_banda_sim (supplementary)
)

for (s in scripts_en) {
  cat(sprintf("\n>> %s\n", s))
  source(s, local = new.env())
}
cat("\n>> done. English figures in", DIR_FIGURAS_EN, "\n")
