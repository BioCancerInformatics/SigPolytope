#################### PART 1 ####################
#' Build geometric latent representation for the main signature
#'
#' @param data Tibble/data.frame com todas as colunas disponíveis.
#' @param cols Vetor de nomes de colunas de `data` que serão usadas
#'   como coordenadas da geometria da assinatura (todas devem ser numéricas
#'   ou coercíveis para numeric).
#' @param id_col Nome da coluna de ID dos circuitries (ex: "Circuitries_id").
#'               Essa coluna NÃO entra na geometria; ela só identifica as linhas.
#' @param meta_cols Vetor de nomes de colunas de metadados a carregar junto
#'   (ex: c("Pathways","Metabolism")). Se NULL, não traz metadados extras
#'   além de `id_col`, se ele existir.
#'
#' @return Lista com:
#'   - features: matriz N x P com as colunas escolhidas
#'   - meta: tibble com os metadados (incluindo id_col, se existir)
#'   - dim_names: nomes das dimensões usadas (colnames(features))
#' @export
build_signature_geometry <- function(data,
                                     cols,
                                     id_col = "Circuitries_id",
                                     meta_cols = NULL) {
  stopifnot(is.data.frame(data))

  # Colunas que realmente existem em `data`
  cols_ok <- intersect(cols, colnames(data))
  if (length(cols_ok) == 0L) {
    stop("Nenhuma coluna de 'cols' foi encontrada em 'data'.")
  }

  # Extrai features e força numeric
  feats <- data[, cols_ok, drop = FALSE]
  feats[] <- lapply(feats, function(x) {
    if (is.numeric(x)) return(x)
    suppressWarnings(as.numeric(x))
  })
  feats <- as.matrix(feats)

  # Metadados
  meta_cols_all <- unique(c(id_col, meta_cols))
  meta_cols_all <- intersect(meta_cols_all, colnames(data))

  if (length(meta_cols_all) > 0L) {
    meta <- tibble::as_tibble(data[, meta_cols_all, drop = FALSE])
  } else {
    meta <- tibble::tibble()
  }

  # rownames = ID (se existir)
  if (id_col %in% colnames(data)) {
    rownames(feats) <- data[[id_col]]
    if (nrow(meta) > 0L) {
      rownames(meta) <- data[[id_col]]
    }
  }

  structure(
    list(
      features  = feats,
      meta      = meta,
      dim_names = colnames(feats)
    ),
    class = "signature_geometry"
  )
}



#################### PART 2 ####################
#' Build geometric latent representation for the regulatory signature
#'
#' Idêntico a build_signature_geometry(), mas semanticamente voltado
#' ao lado regulador (o usuário escolhe outras colunas em `cols`).
#'
#' @inheritParams build_signature_geometry
#'
#' @return Lista análoga à de build_signature_geometry(), mas para o regulador.
#' @export
build_regulator_geometry <- function(data,
                                     cols,
                                     id_col = "Circuitries_id",
                                     meta_cols = NULL) {

  obj <- build_signature_geometry(
    data      = data,
    cols      = cols,
    id_col    = id_col,
    meta_cols = meta_cols
  )

  class(obj) <- c("regulator_geometry", class(obj))
  obj

}


