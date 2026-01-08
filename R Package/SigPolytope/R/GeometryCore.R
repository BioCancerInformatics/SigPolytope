# ===================== GeometryCore.R =====================

#' Internal: PCA embedding from generic tensor (sig + int)
#'
#' @param tensor Lista com:
#'   - features_sig: matriz N x P (assinatura principal)
#'   - features_int: matriz N x P (regulador)
#' @param n_components Número de componentes principais (default = 3).
#'
#' @return Lista com:
#'   - pca: objeto prcomp
#'   - coords_sig: coordenadas N x n_components para o lado sig
#'   - coords_int: coordenadas N x n_components para o lado int
#' @keywords internal


# ===================== GeometryCore.R (PATCH) =====================

#' Internal: project arbitrary matrix into existing PCA
#'
#' @param pca_obj Objeto prcomp com $rotation, $center, $scale.
#' @param X Matriz N x P a ser projetada.
#' @param n_components Número de PCs a retornar.
#' @keywords internal
project_to_pca <- function(pca_obj, X, n_components = 3) {
  X <- as.matrix(X)
  d <- ncol(X)

  if (is.null(pca_obj$rotation)) {
    stop("Objeto PCA sem 'rotation'.")
  }

  # ✅ Correção: rotation é P x K, então a compatibilidade deve checar nrow(rotation) == P
  if (nrow(pca_obj$rotation) != d) {
    stop("Dimensão incompatível: X tem ", d,
         " colunas, rotation tem ", nrow(pca_obj$rotation), " linhas.")
  }

  X_scaled <- scale(
    X,
    center = pca_obj$center,
    scale  = pca_obj$scale
  )

  scores <- X_scaled %*% pca_obj$rotation
  scores[, seq_len(min(n_components, ncol(scores))), drop = FALSE]
}


#' Compute convergence/divergence between signature and regulator geometries
#'
#' @param sig_geom Objeto retornado por build_signature_geometry().
#' @param reg_geom Objeto retornado por build_regulator_geometry().
#'        Deve ter as mesmas linhas (mesmos circuitries, mesma ordem).
#' @param n_components Número de PCs a usar (default 3).
#' @param parallel Logical; se TRUE tenta paralelizar o cálculo dos hulls.
#' @param n_cores Número de cores para paralelismo (default = parallel::detectCores() - 1).
#'
#' @return Uma lista com:
#'   - results: tibble com uma linha por circuitry contendo:
#'       * barycenter_distance
#'       * sig_hull_vol, int_hull_vol
#'       * vol_ratio
#'       * distance_implication, vol_implication
#'       * metadados (Pathways, Metabolism, etc., se estavam em meta)
#'   - embedding: lista com PCA e coords
#'   - pca: objeto prcomp
#' @export
compute_circuitry_convergence <- function(sig_geom,
                                          reg_geom,
                                          n_components = 3,
                                          parallel = FALSE,
                                          n_cores = NULL,
                                          strict_dims = FALSE) {

  stopifnot(is.matrix(sig_geom$features), is.matrix(reg_geom$features))
  stopifnot(nrow(sig_geom$features) == nrow(reg_geom$features))

  # ✅ NEW: guarantee compatible dimensions/order
  aligned <- .align_sig_int_features(sig_geom, reg_geom, strict = strict_dims)

  tensor <- list(
    features_sig = aligned$X_sig,
    features_int = aligned$X_int,
    meta         = sig_geom$meta
  )

  if (!is.null(aligned$dim_names)) {
    tensor$latent_dims <- list(
      canonical = aligned$dim_names,
      sig = aligned$dim_names,
      int = aligned$dim_names
    )
  }

  embedding <- compute_pca_embedding_from_tensor(tensor, n_components = n_components)
  # --------------------- helper para extrair volume ----------------------
  extract_volume <- function(h) {
    # Caso nulo
    if (is.null(h)) return(NA_real_)

    # Se for lista do geometry::convhulln
    if (is.list(h)) {
      if (!is.null(h$volume) && length(h$volume) == 1L) {
        return(as.numeric(h$volume))
      }
      if (!is.null(h$vol) && length(h$vol) == 1L) {
        return(as.numeric(h$vol))
      }
    }

    # Se já for um número
    if (is.numeric(h) && length(h) == 1L) {
      return(as.numeric(h))
    }

    # Qualquer outra coisa → NA
    NA_real_
  }

  # ------------------- helper para 1 circuitry ---------------------------
  compute_geometry_for_index <- function(i) {
    poly_i <- build_circuitry_polytope(tensor, embedding, index = i)

    # distância entre baricentros (sempre deve existir)
    d_bary <- sqrt(sum((poly_i$bary_sig_3d - poly_i$bary_int_3d)^2))

    # tentar extrair volume dos hulls
    v_sig <- extract_volume(poly_i$hull_sig)
    v_int <- extract_volume(poly_i$hull_int)

    # blindagem: se vier numeric(0) ou não finito, vira NA
    if (length(v_sig) == 0L || !is.finite(v_sig)) v_sig <- NA_real_
    if (length(v_int) == 0L || !is.finite(v_int)) v_int <- NA_real_

    list(
      barycenter_distance = d_bary,
      sig_hull_vol        = v_sig,
      int_hull_vol        = v_int
    )
  }
  # ----------------------------------------------------------------------

  # ----------------------------------------------------------------------

  # 4) Calcular distância + volumes para todos os circuitries
  N       <- nrow(tensor$features_sig)
  indices <- seq_len(N)

  if (isTRUE(parallel)) {
    if (is.null(n_cores)) {
      n_cores <- max(1L, parallel::detectCores() - 1L)
    }
    cl <- parallel::makeCluster(n_cores)

    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(geometry)
        library(tibble)
      })
    })

    parallel::clusterExport(
      cl,
      varlist = c("tensor", "embedding", "build_circuitry_polytope", "extract_volume"),
      envir   = environment()
    )

    res_list <- parallel::parLapply(cl, indices, compute_geometry_for_index)
    parallel::stopCluster(cl)
  } else {
    res_list <- lapply(indices, compute_geometry_for_index)
  }

  # 5) Desempacotar
  barycenter_distance_vec <- vapply(res_list, `[[`, numeric(1), "barycenter_distance")
  sig_hull_vol_vec        <- vapply(res_list, `[[`, numeric(1), "sig_hull_vol")
  int_hull_vol_vec        <- vapply(res_list, `[[`, numeric(1), "int_hull_vol")

  V_mean_vec <- rowMeans(cbind(sig_hull_vol_vec, int_hull_vol_vec), na.rm = TRUE)
  eps_vol    <- 1e-12
  vol_min    <- pmin(sig_hull_vol_vec, int_hull_vol_vec)
  vol_max    <- pmax(sig_hull_vol_vec, int_hull_vol_vec)
  vol_ratio_vec <- vol_max / pmax(vol_min, eps_vol)
  vol_ratio_vec[!is.finite(vol_ratio_vec)] <- NA_real_

  # 6) Classificação geométrica
  distance_breaks <- c(0.5, 1.5, 2.5)
  finite_vm <- V_mean_vec[is.finite(V_mean_vec) & V_mean_vec > 0]
  if (length(finite_vm) >= 10L) {
    vol_breaks <- as.numeric(stats::quantile(finite_vm, probs = c(1/3, 2/3)))
  } else {
    vm_med     <- stats::median(finite_vm, na.rm = TRUE)
    vol_breaks <- c(vm_med / 2, vm_med * 2)
  }

  interpret_geometric_pattern <- function(d_bary, vol_ratio, v_mean) {
    distance_implication <- dplyr::case_when(
      !is.finite(d_bary)                  ~ "distance_undefined",
      d_bary < distance_breaks[1]         ~ "high_concordance",
      d_bary < distance_breaks[2]         ~ "moderate_discordance",
      d_bary < distance_breaks[3]         ~ "strong_discordance",
      TRUE                                ~ "extreme_discordance"
    )

    if (!is.finite(v_mean) || v_mean <= 0) {
      vol_implication <- "flat_or_unresolved_geometry"
    } else {
      vol_complexity <- dplyr::case_when(
        v_mean <= vol_breaks[1] ~ "low_dimensional_flat",
        v_mean <= vol_breaks[2] ~ "intermediate_complexity",
        TRUE                    ~ "high_complexity_multidimensional"
      )
      vol_asymmetry <- dplyr::case_when(
        !is.finite(vol_ratio)            ~ "symmetric_or_undefined",
        vol_ratio < 1.5                  ~ "symmetric_sig_int",
        vol_ratio < 3                    ~ "moderately_asymmetric_sig_int",
        TRUE                             ~ "strongly_asymmetric_sig_int"
      )
      vol_implication <- paste(vol_complexity, vol_asymmetry, sep = " ; ")
    }

    list(
      distance_implication = distance_implication,
      vol_implication      = vol_implication
    )
  }

  imp_list <- Map(
    interpret_geometric_pattern,
    d_bary    = barycenter_distance_vec,
    vol_ratio = vol_ratio_vec,
    v_mean    = V_mean_vec
  )

  distance_implication_vec <- vapply(imp_list, `[[`, character(1), "distance_implication")
  vol_implication_vec      <- vapply(imp_list, `[[`, character(1), "vol_implication")

  # 7) Montar tibble final com metadados
  results <- dplyr::bind_cols(
    sig_geom$meta,
    tibble::tibble(
      barycenter_distance  = barycenter_distance_vec,
      distance_implication = distance_implication_vec,
      sig_hull_vol         = sig_hull_vol_vec,
      int_hull_vol         = int_hull_vol_vec,
      vol_ratio            = vol_ratio_vec,
      vol_implication      = vol_implication_vec
    )
  )

  list(
    results   = results,
    embedding = embedding,
    pca       = embedding$pca
  )
}



#' Internal: build barycenter-centered +/- vertices for one latent vector
#'
#' @param latent_vec Vetor numérico de dimensão P.
#' @param side Texto "sig" ou "int".
#' @param pca_obj Objeto prcomp (da compute_pca_embedding_from_tensor).
#' @param n_components Número de PCs.
#' @param eps Passo mínimo ao redor de zero.
#' @keywords internal
build_vertices_for_signature <- function(latent_vec,
                                         side = c("sig", "int"),
                                         pca_obj,
                                         n_components = 3,
                                         eps = 1e-6) {
  side <- match.arg(side)
  v <- as.numeric(latent_vec)
  P <- length(v)

  step <- pmax(abs(v), eps)

  verts_latent <- matrix(0, nrow = 2 * P, ncol = P)
  dim_idx      <- integer(2 * P)
  sign_flag    <- character(2 * P)

  row_idx <- 1L
  for (j in seq_len(P)) {
    # +step
    vp <- v
    vp[j] <- v[j] + step[j]
    verts_latent[row_idx, ] <- vp
    dim_idx[row_idx]        <- j
    sign_flag[row_idx]      <- "pos"
    row_idx <- row_idx + 1L

    # -step
    vn <- v
    vn[j] <- v[j] - step[j]
    verts_latent[row_idx, ] <- vn
    dim_idx[row_idx]        <- j
    sign_flag[row_idx]      <- "neg"
    row_idx <- row_idx + 1L
  }

  verts_pca <- project_to_pca(pca_obj, verts_latent, n_components = n_components)
  colnames(verts_pca) <- paste0("PC", seq_len(ncol(verts_pca)))

  tibble::tibble(
    side      = side,
    dim_index = dim_idx,
    sign      = sign_flag,
    PC1       = verts_pca[, 1],
    PC2       = verts_pca[, 2],
    PC3       = verts_pca[, 3]
  )
}

#' Internal: build 3D polytope for one circuitry (sig + int)
#'
#' @param tensor Lista com features_sig, features_int, meta.
#' @param embedding Lista retornada por compute_pca_embedding_from_tensor().
#' @param index Índice da linha (circuitry) a ser usada.
#'
#' @return Lista com schema canônico esperado pelos plots:
#'   - pts_sig_3d / pts_int_3d (matrizes 3D)
#'   - bary_sig_3d / bary_int_3d
#'   - hull_sig / hull_int
#'   - (opcional) vertices_* para debug
#' @keywords internal
build_circuitry_polytope <- function(tensor,
                                     embedding,
                                     index) {

  meta <- tensor$meta
  n_meta <- if (!is.null(meta)) nrow(meta) else nrow(tensor$features_sig)
  if (index < 1L || index > n_meta) stop("Index fora do intervalo.")

  v_sig <- tensor$features_sig[index, ]
  v_int <- tensor$features_int[index, ]

  pca_obj <- embedding$pca
  if (is.null(pca_obj)) {
    stop("embedding$pca é NULL.")
  }

  # nomes das dimensões (genérico, sem 18D fixo)
  dim_names <- colnames(tensor$features_sig)

  # vértices (em PCA)
  verts_sig <- build_vertices_for_signature(
    latent_vec   = v_sig,
    side         = "sig",
    pca_obj      = pca_obj,
    n_components = 3
  )
  verts_int <- build_vertices_for_signature(
    latent_vec   = v_int,
    side         = "int",
    pca_obj      = pca_obj,
    n_components = 3
  )

  verts_sig$dim_name <- dim_names[verts_sig$dim_index]
  verts_int$dim_name <- dim_names[verts_int$dim_index]

  vertices_all <- dplyr::bind_rows(verts_sig, verts_int)

  # ✅ MATRIZES 3D CANÔNICAS ESPERADAS PELOS PLOTS
  pts_sig_3d <- as.matrix(verts_sig[, c("PC1", "PC2", "PC3")])
  pts_int_3d <- as.matrix(verts_int[, c("PC1", "PC2", "PC3")])

  # barycenters (projeção 3D)
  bary_sig_3d <- project_to_pca(pca_obj, matrix(v_sig, nrow = 1), n_components = 3)[1, ]
  bary_int_3d <- project_to_pca(pca_obj, matrix(v_int, nrow = 1), n_components = 3)[1, ]
  names(bary_sig_3d) <- c("PC1", "PC2", "PC3")
  names(bary_int_3d) <- c("PC1", "PC2", "PC3")

  # convex hull sig
  hull_sig <- list(faces = NULL, volume = NA_real_, area = NA_real_)
  if (nrow(pts_sig_3d) >= 4L) {
    Hs <- try(geometry::convhulln(pts_sig_3d, options = "FA"), silent = TRUE)
    if (!inherits(Hs, "try-error")) {
      faces_s <- Hs
      vol_s   <- attr(Hs, "vol")
      area_s  <- attr(Hs, "area")
      if (is.list(Hs)) {
        faces_s <- if (!is.null(Hs$hull)) Hs$hull else if (!is.null(Hs$facets)) Hs$facets else NULL
        if (!is.null(Hs$vol))  vol_s  <- Hs$vol
        if (!is.null(Hs$area)) area_s <- Hs$area
      }
      hull_sig$faces  <- faces_s
      hull_sig$volume <- vol_s
      hull_sig$area   <- area_s
    }
  }

  # convex hull int
  hull_int <- list(faces = NULL, volume = NA_real_, area = NA_real_)
  if (nrow(pts_int_3d) >= 4L) {
    Hi <- try(geometry::convhulln(pts_int_3d, options = "FA"), silent = TRUE)
    if (!inherits(Hi, "try-error")) {
      faces_i <- Hi
      vol_i   <- attr(Hi, "vol")
      area_i  <- attr(Hi, "area")
      if (is.list(Hi)) {
        faces_i <- if (!is.null(Hi$hull)) Hi$hull else if (!is.null(Hi$facets)) Hi$facets else NULL
        if (!is.null(Hi$vol))  vol_i  <- Hi$vol
        if (!is.null(Hi$area)) area_i <- Hi$area
      }
      hull_int$faces  <- faces_i
      hull_int$volume <- vol_i
      hull_int$area   <- area_i
    }
  }

  circuitry_id_out <- if (!is.null(meta) && "Circuitries_id" %in% colnames(meta)) {
    meta$Circuitries_id[index]
  } else {
    as.character(index)
  }

  # ✅ Retorno canônico para os plots (mantendo vertices_* como debug)
  list(
    circuitry_id = circuitry_id_out,

    pts_sig_3d   = pts_sig_3d,
    pts_int_3d   = pts_int_3d,

    bary_sig_3d  = bary_sig_3d,
    bary_int_3d  = bary_int_3d,

    hull_sig     = hull_sig,
    hull_int     = hull_int,

    # (opcional) debug/inspeção
    vertices_sig = verts_sig,
    vertices_int = verts_int,
    vertices_all = vertices_all
  )
}
