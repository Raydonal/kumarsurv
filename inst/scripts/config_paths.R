## config_paths.R
## Resolucao de caminhos do PACOTE_VALIDACAO.
##
## MOTIVO (Rodada 33): os scripts usavam tres convencoes incompativeis --
## dados/arq.csv, ../dados/arq.csv e nomes nus para os .rds, que na verdade
## vivem em resultados/ -- alem de gravar figuras no caminho absoluto
## /home/claude/build/Imgs/, inexistente fora do ambiente original. Nenhum
## script rodava a partir de scripts/ sem preparacao manual.
##
## SOLUCAO: localizar a raiz do pacote e instalar resolvedores que aceitam
## qualquer uma das convencoes antigas. Leitura: se o caminho dado existir, e'
## usado sem alteracao (nada quebra); senao, procura-se o nome-base no
## diretorio canonico. Escrita: se o diretorio do caminho dado nao existir,
## grava-se no diretorio canonico. A logica numerica dos scripts nao e' tocada.
##
## Sobrescrever a raiz, se necessario:  Sys.setenv(PACOTE_VALIDACAO = "/caminho")

if (!exists(".PATHS_CONFIGURADOS", envir = globalenv())) {

  .achar_raiz <- function() {
    env <- Sys.getenv("PACOTE_VALIDACAO", "")
    if (nzchar(env) && dir.exists(env)) return(normalizePath(env))
    ## sobe a partir do diretorio de trabalho procurando a marca do pacote
    d <- normalizePath(getwd())
    for (i in 1:6) {
      if (dir.exists(file.path(d, "dados")) &&
          dir.exists(file.path(d, "scripts"))) return(d)
      pai <- dirname(d)
      if (identical(pai, d)) break
      d <- pai
    }
    ## fallback: diretorio de trabalho corrente (modo "tudo lado a lado")
    normalizePath(getwd())
  }

  RAIZ           <- .achar_raiz()
  DIR_DADOS      <- if (dir.exists(file.path(RAIZ, "dados")))      file.path(RAIZ, "dados")      else RAIZ
  DIR_RESULTADOS <- if (dir.exists(file.path(RAIZ, "resultados"))) file.path(RAIZ, "resultados") else RAIZ
  DIR_SCRIPTS    <- if (dir.exists(file.path(RAIZ, "scripts")))    file.path(RAIZ, "scripts")    else RAIZ
  DIR_FIGURAS    <- file.path(RAIZ, "figuras")
  if (!dir.exists(DIR_FIGURAS)) dir.create(DIR_FIGURAS, recursive = TRUE, showWarnings = FALSE)

  ## --- resolvedores -------------------------------------------------------
  ## leitura: preserva o caminho se ele existir; senao procura o nome-base
  ## no diretorio canonico e, em ultimo caso, nos demais diretorios do pacote.
  .resolver_leitura <- function(f, dir_canonico) {
    if (file.exists(f)) return(f)
    cand <- file.path(dir_canonico, basename(f))
    if (file.exists(cand)) return(cand)
    for (d in c(DIR_RESULTADOS, DIR_DADOS, DIR_SCRIPTS, RAIZ)) {
      cand <- file.path(d, basename(f))
      if (file.exists(cand)) return(cand)
    }
    f  # deixa falhar com a mensagem original do R
  }

  ## escrita: se o diretorio pedido existir, respeita-o; senao usa o canonico.
  .resolver_escrita <- function(f, dir_canonico) {
    d <- dirname(f)
    if (d != "." && dir.exists(d)) return(f)
    file.path(dir_canonico, basename(f))
  }

  ## --- shims ---------------------------------------------------------------
  readRDS <- function(file, ...) base::readRDS(.resolver_leitura(file, DIR_RESULTADOS), ...)
  saveRDS <- function(object, file, ...) base::saveRDS(object, .resolver_escrita(file, DIR_RESULTADOS), ...)
  read.csv <- function(file, ...) utils::read.csv(.resolver_leitura(file, DIR_DADOS), ...)
  read.table <- function(file, ...) utils::read.table(.resolver_leitura(file, DIR_DADOS), ...)
  write.csv <- function(x, file = "", ...) utils::write.csv(x, .resolver_escrita(file, DIR_RESULTADOS), ...)

  cairo_pdf <- function(filename = "Rplot.pdf", ...)
    grDevices::cairo_pdf(.resolver_escrita(filename, DIR_FIGURAS), ...)
  pdf <- function(file = "Rplot.pdf", ...)
    grDevices::pdf(.resolver_escrita(file, DIR_FIGURAS), ...)
  png <- function(filename = "Rplot.png", ...)
    grDevices::png(.resolver_escrita(filename, DIR_FIGURAS), ...)

  ## ggsave so' e' redefinido se o ggplot2 estiver disponivel
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggsave <- function(filename, ...)
      ggplot2::ggsave(.resolver_escrita(filename, DIR_FIGURAS), ...)
  }

  .PATHS_CONFIGURADOS <- TRUE
  if (!identical(Sys.getenv("PACOTE_VALIDACAO_QUIET"), "1")) {
    cat(sprintf("[config_paths] raiz=%s | dados=%s | resultados=%s | figuras=%s\n",
                RAIZ, basename(DIR_DADOS), basename(DIR_RESULTADOS), basename(DIR_FIGURAS)))
  }
}
