#' Build geometry from (direction, strength) pairs
#'
#' @param data data.frame/tibble com os dados.
#' @param dir_cols Named character vector:
#'        names = nomes das dimensões latentes (ex: "rho", "TN", "OS")
#'        values = nomes das colunas de direção em `data`
#' @param strength_cols Named character vector com MESMOS names de dir_cols:
#'        values = nomes das colunas de força em `data`
#' @param id_col Coluna de ID
#' @param meta_cols Metadados
#' @return Objeto da classe `signature_geometry`
#' @export
build_signature_geometry_from_pairs <- function(
    data,
    dir_cols,
    strength_cols,
    id_col = "Circuitries_id",
    meta_cols = NULL
) {
  stopifnot(is.data.frame(data))

  # Garante alinhamento de nomes
  if (!setequal(names(dir_cols), names(strength_cols))) {
    stop("dir_cols e strength_cols devem ter os MESMOS nomes (dimensões latentes).")
  }

  if (length(dir_cols) < 3L) {
    stop("Você precisa de pelo menos 3 pares (direção + força) para construir a geometria.")
  }

  dims <- names(dir_cols)

  encode_dir_strength <- function(dir_vec, strength_vec) {
    f <- factor(dir_vec)
    lv <- levels(f)
    k  <- length(lv)

    if (k == 1L) {
      dir_code <- rep(0, length(dir_vec))
    } else if (k == 2L) {
      dir_code <- ifelse(f == lv[1], -1, ifelse(f == lv[2], 1, 0))
    } else {
      dir_code <- seq(-1, 1, length.out = k)[as.integer(f)]
    }

    as.numeric(dir_code) * suppressWarnings(as.numeric(strength_vec))
  }

  feat_list <- lapply(dims, function(dim_name) {
    dir_col <- dir_cols[[dim_name]]
    str_col <- strength_cols[[dim_name]]

    if (!dir_col %in% colnames(data)) stop("Coluna de direção não encontrada: ", dir_col)
    if (!str_col %in% colnames(data)) stop("Coluna de força não encontrada: ", str_col)

    encode_dir_strength(
      dir_vec      = data[[dir_col]],
      strength_vec = data[[str_col]]
    )
  })

  feats <- as.data.frame(feat_list, optional = TRUE, stringsAsFactors = FALSE)
  colnames(feats) <- dims
  feats <- as.matrix(feats)

  # Metadados
  meta_cols_all <- unique(c(id_col, meta_cols))
  meta_cols_all <- intersect(meta_cols_all, colnames(data))

  if (length(meta_cols_all) > 0L) {
    meta <- tibble::as_tibble(data[, meta_cols_all, drop = FALSE])
  } else {
    meta <- tibble::tibble()
  }

  # rownames = ID
  if (id_col %in% colnames(data)) {
    rownames(feats) <- data[[id_col]]
    if (nrow(meta) > 0L) {
      meta_df <- as.data.frame(meta)
      rownames(meta_df) <- data[[id_col]]
      meta <- tibble::as_tibble(meta_df)
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

#' Build geometry from (direction, strength) pairs for the regulator side
#'
#' @inheritParams build_signature_geometry_from_pairs
#' @return Objeto da classe `regulator_geometry`
#' @export
build_regulator_geometry_from_pairs <- function(
    data,
    dir_cols,
    strength_cols,
    id_col = "Circuitries_id",
    meta_cols = NULL
) {
  obj <- build_signature_geometry_from_pairs(
    data          = data,
    dir_cols      = dir_cols,
    strength_cols = strength_cols,
    id_col        = id_col,
    meta_cols     = meta_cols
  )

  class(obj) <- c("regulator_geometry", class(obj))
  obj
}
