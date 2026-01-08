# ============================================================
# R/align_features.R
# ============================================================
# Internal: guarantee sig/int feature matrices are compatible
# (same dims, same order). Global: relies on column names when available.
# ============================================================

#' @keywords internal
.align_sig_int_features <- function(sig_geom, reg_geom, strict = TRUE) {
  stopifnot(is.matrix(sig_geom$features), is.matrix(reg_geom$features))
  
  Xs <- sig_geom$features
  Xi <- reg_geom$features
  
  # 0) same number of circuitries (rows)
  if (nrow(Xs) != nrow(Xi)) {
    stop("sig_geom e reg_geom devem ter o mesmo número de linhas (mesmos circuitries).")
  }
  
  ns <- colnames(Xs)
  ni <- colnames(Xi)
  
  # 1) If both have column names -> align by intersection + identical order
  if (!is.null(ns) && !is.null(ni)) {
    
    common <- intersect(ns, ni)
    
    if (length(common) == 0L) {
      msg <- paste0(
        "sig_geom$features e reg_geom$features não possuem dimensões em comum (por colnames).\n",
        "Dica: garanta colnames consistentes (ex.: 'TN_dir', 'OS_dir', etc.) ou use build_tensor_semantic."
      )
      stop(msg)
    }
    
    # Keep sig order as canonical (you can flip to schema order if you prefer)
    Xs2 <- Xs[, common, drop = FALSE]
    Xi2 <- Xi[, common, drop = FALSE]
    
    # Ensure same order
    Xi2 <- Xi2[, colnames(Xs2), drop = FALSE]
    
    # Optional strictness: require no dimension loss
    if (strict) {
      if (length(common) != length(ns) || length(common) != length(ni)) {
        missing_sig <- setdiff(ns, common)
        missing_int <- setdiff(ni, common)
        stop(
          "Dimensões incompatíveis entre sig/int.\n",
          "Faltando em reg (presentes em sig): ", paste(missing_sig, collapse = ", "), "\n",
          "Faltando em sig (presentes em reg): ", paste(missing_int, collapse = ", ")
        )
      }
    }
    
    return(list(
      X_sig = Xs2,
      X_int = Xi2,
      dim_names = colnames(Xs2)
    ))
  }
  
  # 2) Without column names -> we can only check shape
  if (ncol(Xs) != ncol(Xi)) {
    stop(
      "Dimensões incompatíveis: sig possui ", ncol(Xs),
      " colunas e reg possui ", ncol(Xi), " colunas.\n",
      "Dica: adicione colnames ou use build_tensor_semantic para padronizar dimensões."
    )
  }
  
  # assume already aligned
  list(
    X_sig = Xs,
    X_int = Xi,
    dim_names = NULL
  )
}
