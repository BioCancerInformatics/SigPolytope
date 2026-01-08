# ============================================================
# PCA utils for GeoCircuitry
# ============================================================

#' PCA embedding from two feature matrices (signature x regulator)
#'
#' @param tensor Lista com:
#'   - features_sig: matriz N x P
#'   - features_int: matriz N x P
#' @param n_components Número de PCs a manter (default = 3)
#'
#' @return Lista com:
#'   - pca: objeto prcomp
#'   - coords_sig: N x n_components (scores da assinatura)
#'   - coords_int: N x n_components (scores do regulador)
#' @keywords internal
compute_pca_embedding_from_tensor <- function(tensor, n_components = 3) {
  stopifnot(
    !is.null(tensor$features_sig),
    !is.null(tensor$features_int)
  )

  X_sig <- tensor$features_sig
  X_int <- tensor$features_int

  if (!is.matrix(X_sig) || !is.matrix(X_int)) {
    stop("tensor$features_sig e tensor$features_int devem ser matrizes numéricas.")
  }
  if (ncol(X_sig) != ncol(X_int)) {
    stop("Sig e Int precisam ter o mesmo número de colunas (dimensões latentes).")
  }

  # Empilha as duas metades (2N x P)
  X_all <- rbind(X_sig, X_int)

  # Estatísticas de escala
  means <- colMeans(X_all, na.rm = TRUE)
  sds   <- apply(X_all, 2, sd, na.rm = TRUE)

  # Corrige sd = 0 ou NA
  zero_var_idx <- which(!is.finite(sds) | sds == 0)
  if (length(zero_var_idx) > 0L) {
    warning(
      "Dimensões sem variância ou com sd não finito: ",
      paste(colnames(X_all)[zero_var_idx], collapse = ", "),
      ". Ajuste: sd = 1."
    )
    sds[zero_var_idx] <- 1
  }

  # Padronização
  X_scaled <- sweep(X_all, 2, means, "-")
  X_scaled <- sweep(X_scaled, 2, sds,   "/")

  # PCA global
  pca_obj <- stats::prcomp(X_scaled, center = FALSE, scale. = FALSE)
  pca_obj$center <- means
  pca_obj$scale  <- sds

  # Scores
  k <- min(n_components, ncol(pca_obj$x))
  coords_all <- pca_obj$x[, seq_len(k), drop = FALSE]

  N <- nrow(X_sig)
  coords_sig <- coords_all[1:N, , drop = FALSE]
  coords_int <- coords_all[(N + 1):(2 * N), , drop = FALSE]

  colnames(coords_sig) <- paste0("PC", seq_len(ncol(coords_sig)))
  colnames(coords_int) <- paste0("PC", seq_len(ncol(coords_int)))

  list(
    pca        = pca_obj,
    coords_sig = coords_sig,
    coords_int = coords_int
  )
}



#' Project new points into an existing PCA space
#'
#' @param pca_obj Objeto prcomp retornado por compute_pca_embedding_from_tensor()
#' @param X Matriz N x P com as mesmas colunas usadas no PCA
#' @param n_components Número de PCs a retornar
#'
#' @return Matriz N x n_components com novas coordenadas projetadas
#' @keywords internal
project_to_pca <- function(pca_obj, X, n_components = 3) {
  X <- as.matrix(X)
  d <- ncol(X)

  if (is.null(pca_obj$rotation)) {
    stop("Objeto PCA não contém 'rotation'. Verifique o PCA gerado.")
  }

  # A CHECAGEM CORRETA É COM nrow(rotation)
  if (nrow(pca_obj$rotation) != d) {
    stop(
      "Dimensão incompatível: X tem ", d,
      " colunas, mas rotation tem ", nrow(pca_obj$rotation), " linhas."
    )
  }

  # Padronização consistente
  X_scaled <- scale(
    X,
    center = pca_obj$center,
    scale  = pca_obj$scale
  )

  # Quantos PCs projetar
  k <- min(n_components, ncol(pca_obj$rotation))

  scores <- X_scaled %*% pca_obj$rotation[, seq_len(k), drop = FALSE]
  colnames(scores) <- paste0("PC", seq_len(k))

  scores
}
