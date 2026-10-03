## run_all.R — reproduz a analise do artigo, na ordem de dependencia.
##
## Rodada 33: nao exige mais montar um diretorio unico "tudo lado a lado".
## config_paths.R localiza a raiz do pacote e resolve dados/, resultados/ e
## figuras/ a partir de qualquer diretorio de trabalho. Para rodar:
##
##   LC_ALL=C.utf8 LANG=C.utf8 Rscript scripts/run_all.R
##   LC_ALL=C.utf8 LANG=C.utf8 Rscript run_all.R          # de dentro de scripts/
##
## Raiz explicita, se necessario:
##   PACOTE_VALIDACAO=/caminho/para/PACOTE_VALIDACAO Rscript scripts/run_all.R
##
## Selecao dos blocos por variavel de ambiente (padrao: rapidos + figuras):
##   RUN_PESADOS=1   inclui os experimentos com M=300 / B=300 (cerca de 3 h no total,
##                   medido na reexecucao de 14-15/09/2026 com a biblioteca da Rodada 47)
##   RUN_ONLY=03,06  roda apenas os scripts cujos prefixos forem listados

## Rodada 41: o source() abaixo era relativo ao diretorio corrente, de modo que
## a forma documentada "Rscript scripts/run_all.R" a partir da raiz do pacote
## falhava antes de qualquer coisa. Localiza-se scripts/ primeiro.
if (!file.exists("copula_ml.R")) {
  .p <- Sys.getenv("PACOTE_VALIDACAO", "")
  if (!nzchar(.p) || !dir.exists(.p)) {
    .p <- normalizePath(getwd())
    for (.i in 1:6) {
      if (dir.exists(file.path(.p, "scripts")) && dir.exists(file.path(.p, "dados"))) break
      .q <- dirname(.p); if (identical(.q, .p)) break; .p <- .q
    }
  }
  if (dir.exists(file.path(.p, "scripts"))) setwd(file.path(.p, "scripts"))
}
source("copula_ml.R")   # carrega config_paths.R
setwd(DIR_SCRIPTS)      # todos os source() abaixo usam nomes nus

.only <- Sys.getenv("RUN_ONLY", "")
.only <- if (nzchar(.only)) trimws(strsplit(.only, ",")[[1]]) else character(0)
.pesados <- identical(Sys.getenv("RUN_PESADOS"), "1")

.rodar <- function(arquivo, descricao, pesado = FALSE) {
  if (length(.only) && !any(startsWith(arquivo, .only))) return(invisible(NULL))
  if (pesado && !.pesados) {
    cat(sprintf(">> [pulado, defina RUN_PESADOS=1] %s\n", descricao)); return(invisible(NULL))
  }
  cat(sprintf("\n>> %s  (%s)\n", descricao, arquivo))
  t0 <- Sys.time()
  ok <- tryCatch({ source(arquivo, local = new.env()); TRUE },
                 error = function(e) { cat("   ERRO: ", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("   %s em %.1f min\n", if (ok) "concluido" else "FALHOU",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  invisible(ok)
}

## --- 1. ajuste nos dados reais: produz reais_result_sem.rds, do qual --------
##        praticamente todo o resto depende ---------------------------------
.rodar("01_dados_reais_SEM.R", "ajustes nos dados reais (Kuma vs beta)")
.rodar("02_farrington_SEM.R",  "Farrington e CUSUM")

## --- 2. analises que consomem o ajuste -------------------------------------
.rodar("03_quantil_SEM.R",            "regressao quantilica (perda pinball)")
.rodar("08_vantagem_SEM.R",           "avaliacao preditiva fora da amostra")
.rodar("tab_arma_selecao.R",          "selecao da ordem ARMA (tab:arma)")
.rodar("perfil_ar1.R",                "intervalo de perfil do AR(1)")
.rodar("perfil_ar1_exato.R",          "perfil exato do AR(1) nas duas bacias e comparacao de periodicidade por modo (Rodada 47)")
.rodar("analise_multiquantil.R",      "tabela multi-quantil (tab:multitau)")
.rodar("calib_SEM.R",                 "calibracao dos limites em toda a Fase II (tab:calib ate o v11)")
.rodar("calib_controle_SEM.R",        "calibracao sob controle e deteccao no surto (Tabela 5 do artigo; Tabelas S12 e S13)")
.rodar("multitau_controle.R",         "tab:multitau avaliada sob controle (Rodada 47)")
.rodar("analise_sensibilidade_zeros.R", "sensibilidade a epsilon (tab:sensib)")
.rodar("sensib_zeros_controle.R",     "sensibilidade a epsilon, calibracao sob controle (tab:sensib; Rodada 47)")
.rodar("analise_harmonicos_inflacao.R", "periodicidade e inflacao de zeros")
.rodar("analise_inflacao_diagnostico.R", "diagnostico residual: Ljung-Box e ACF (Secao 6.2)")
.rodar("analise_covariaveis.R",   "significancia dos coeficientes e especificacoes de tendencia")
.rodar("teste_multi_tau.R",        "validacao do criterio multi-tau (nota da tab:multitau)")

## --- 3. calculos de simulacao leves ----------------------------------------
.rodar("09_efeito_correlacao.R", "efeito da correlacao (analitico, tab:corr)")
.rodar("10_banda_simulada.R",    "cobertura da banda de predicao")
.rodar("11_arl.R",               "calibracao do ARL (fig:arl)")

## --- 4. figuras -------------------------------------------------------------
.rodar("fig_densidades.R",       "densidades e quantis Kumaraswamy")
.rodar("fig3_horizontal.R",      "serie observada")
.rodar("fig4_cusum_SEM.R",       "carta CUSUM")
.rodar("06_figuras_SEM.R",       "carta quantilica: limites marginais (Figura 3 do artigo)")
.rodar("fig_quantis_corrigida.R", "quantis marginais em escala log (Figura S4 do suplementar)")
.rodar("07_fig_deteccao_SEM.R",  "deteccao: CUSUM vs Farrington")

## --- 5. experimentos caros (M=300 / B=300) ----------------------------------
.rodar("13_bootstrap_gof.R",          "bootstrap de bondade de ajuste (B=300)", pesado = TRUE)
.rodar("14_mc_validacao_B300.R",      "Monte Carlo do estimador (M=300)",       pesado = TRUE)
.rodar("16_simcomp_recovery_B300.R",  "recuperacao comparada Kuma/beta (M=300)", pesado = TRUE)
.rodar("17_crps_fix_B300.R",          "custo de ma-especificacao, CRPS (M=300)", pesado = TRUE)
.rodar("15_sim_comp_B300.R",          "perda pinball: direto em tau vs derivado (M=300)", pesado = TRUE)
## Nota (Rodada 41): 15 gera a comparacao PINBALL de sec:sim-comp; seu CRPS foi
## superado por 17_crps_fix_B300.R, que e' a origem dos valores de CRPS no artigo.
.rodar("18_calib_cauda.R",            "calibracao em cauda sob inflacao (M=300)", pesado = TRUE)
.rodar("simulacao_desenho_aplicado_300.R", "recuperacao no desenho da aplicacao (M=300, Tabela S4 do suplementar)", pesado = TRUE)
## Nota (Rodada 41): a transcricao arquivada deste experimento
## (resultados/resultado_simulacao_aplicada_300_FINAL.txt, Quadro 4 do suplementar)
## foi produzida pela via incremental simulacao_checkpoint.R, que acumula
## replicas em sim_checkpoint.rds e permite retomar entre chamadas. O script
## acima faz o mesmo experimento numa unica execucao.
## Rodada 47: a partir desta rodada o Quadro 4 (Tabela S4 desde o supplementary_v13) vem desta execucao unica, com a biblioteca
## corrigida (a transcricao FINAL foi regerada); o experimento abaixo mostra que o AR junto
## a fronteira no Quadro 4 vem dos empates criados pelo truncamento em 1e-7.
.rodar("simaplic_empates.R",          "desenho da aplicacao com e sem truncamento em 1e-7 (Secao S5.1 do suplementar)", pesado = TRUE)

cat(sprintf("\n>> concluido. Figuras em %s, resultados em %s\n", DIR_FIGURAS, DIR_RESULTADOS))
