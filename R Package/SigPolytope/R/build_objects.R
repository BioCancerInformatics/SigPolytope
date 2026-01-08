#' Build a circuitry polytope object (paper-aligned)
#'
#' Separa a etapa de construção (tensor + PCA global + politopo)
#' da etapa de visualização (plot).
#'
#' @param sig_geom Objeto retornado por build_signature_geometry().
#' @param reg_geom Objeto retornado por build_regulator_geometry().
#' @param index Índice do circuitry (1..N).
#' @param n_components Número de PCs a usar (default = 3).
#' @param strict_dims Se TRUE, exige compatibilidade total de dimensões entre sig/int
#'   (sem perda por interseção). Em modo global, o default FALSE é mais tolerante.
#'
#' @return Um objeto (lista) com: tensor, embedding e poly.
#' @export
build_circuitry_polytope_object <- function(sig_geom,
                                            reg_geom,
                                            index = 1,
                                            n_components = 3,
                                            strict_dims = FALSE) {
  stopifnot(
    inherits(sig_geom, "signature_geometry"),
    inherits(reg_geom, "regulator_geometry")
  )
  stopifnot(is.matrix(sig_geom$features), is.matrix(reg_geom$features))

  # --- valida index ---
  if (!is.numeric(index) || length(index) != 1L || is.na(index)) {
    stop("`index` deve ser um inteiro (1..N).")
  }
  index <- as.integer(index)
  if (index < 1L || index > nrow(sig_geom$features)) {
    stop("`index` fora do intervalo: 1..", nrow(sig_geom$features), ".")
  }

  # --- PCA components ---
  if (!is.numeric(n_components) || length(n_components) != 1L || is.na(n_components)) {
    stop("`n_components` deve ser um inteiro >= 2.")
  }
  n_components <- as.integer(n_components)
  if (n_components < 2L) stop("`n_components` deve ser >= 2.")

  # ✅ align dimensions + order (global safety)
  aligned <- .align_sig_int_features(sig_geom, reg_geom, strict = strict_dims)

  tensor <- list(
    features_sig = aligned$X_sig,
    features_int = aligned$X_int,
    meta         = sig_geom$meta
  )

  # attach dims if available (helps annotations downstream)
  if (!is.null(aligned$dim_names)) {
    tensor$latent_dims <- list(
      canonical = aligned$dim_names,
      sig       = aligned$dim_names,
      int       = aligned$dim_names
    )
  }

  embedding <- compute_pca_embedding_from_tensor(
    tensor,
    n_components = n_components
  )

  poly <- build_circuitry_polytope(
    tensor    = tensor,
    embedding = embedding,
    index     = index
  )

  structure(
    list(
      tensor       = tensor,
      embedding    = embedding,
      poly         = poly,
      index        = index,
      n_components = n_components
    ),
    class = "geocircuitry_polytope"
  )
}
