# ============================================================
# Internal preview polytope + safe convex hull (3D)
# ============================================================
# NOTE:
# - build_preview_octahedron() is ONLY a quick preview helper.
# - It does NOT implement the paper method (2P axis-aligned vertices in latent space).
# - .safe_convhull3d() makes hull computation robust to numerical degeneracy,
#   duplicates, NA/Inf, and low-rank (flat) point clouds.
# ============================================================

#' @keywords internal
.safe_convhull3d <- function(coords, joggle = TRUE, tol = 1e-12) {
  coords <- as.matrix(coords)
  if (!is.numeric(coords)) storage.mode(coords) <- "double"

  # keep only finite rows
  coords <- coords[stats::complete.cases(coords), , drop = FALSE]
  if (nrow(coords) == 0) {
    return(list(faces = NULL, volume = 0, area = 0))
  }

  # remove duplicates
  coords <- unique(coords)

  # need at least 4 unique points for a 3D hull
  if (nrow(coords) < 4) {
    return(list(faces = NULL, volume = 0, area = 0))
  }

  # check affine dimension / numerical rank
  X <- scale(coords, center = TRUE, scale = FALSE)
  s <- tryCatch(svd(X, nu = 0, nv = 0)$d, error = function(e) NULL)
  if (is.null(s) || length(s) == 0) {
    return(list(faces = NULL, volume = 0, area = 0))
  }

  # rank threshold relative to max singular value
  thr <- max(s) * tol
  rnk <- sum(s > thr)

  # if rank < 3 => flat/line-like => no 3D volume hull
  if (rnk < 3) {
    return(list(faces = NULL, volume = 0, area = 0))
  }

  # Qhull options
  opt <- if (isTRUE(joggle)) "FA QJ" else "FA"

  H <- try(geometry::convhulln(coords, options = opt), silent = TRUE)

  # fallback: try with QJ if not used
  if (inherits(H, "try-error") && !isTRUE(joggle)) {
    H2 <- try(geometry::convhulln(coords, options = "FA QJ"), silent = TRUE)
    if (!inherits(H2, "try-error")) H <- H2
  }

  if (inherits(H, "try-error")) {
    return(list(faces = NULL, volume = 0, area = 0))
  }

  # geometry::convhulln can return a matrix with attrs (typical)
  # or a list (depending on version/options)
  faces <- NULL
  vol   <- NA_real_
  area  <- NA_real_

  if (is.matrix(H)) {
    faces <- H
    vol   <- suppressWarnings(as.numeric(attr(H, "vol")))
    area  <- suppressWarnings(as.numeric(attr(H, "area")))
  } else if (is.list(H)) {
    faces <- if (!is.null(H$hull)) H$hull else if (!is.null(H$facets)) H$facets else NULL
    vol   <- suppressWarnings(as.numeric(if (!is.null(H$vol)) H$vol else attr(H, "vol")))
    area  <- suppressWarnings(as.numeric(if (!is.null(H$area)) H$area else attr(H, "area")))
  }

  if (!is.finite(vol))  vol  <- 0
  if (!is.finite(area)) area <- 0

  list(faces = faces, volume = vol, area = area)
}

#' Internal: build preview octahedron polytope (NOT paper method)
#'
#' Esta função cria um octaedro 3D simples ao redor dos baricentros,
#' com escala baseada na variabilidade das features. Ela existe apenas
#' como *preview rápido* e NÃO implementa o método do paper (SigPolytope),
#' que usa 2P vértices axis-aligned no espaço latente e convex hull após projeção.
#'
#' @param tensor Lista com features_sig, features_int, meta
#' @param embedding Lista com coords_sig, coords_int
#' @param index Índice do circuitry (1..N)
#' @param scale_factor Fator de escala (multiplica o SD das features do circuitry)
#' @param joggle Se TRUE, usa Qhull "QJ" para evitar falhas numéricas
#'
#' @return Lista simplificada com vértices e hulls (preview)
#' @keywords internal
build_preview_octahedron <- function(tensor,
                                     embedding,
                                     index,
                                     scale_factor = 0.5,
                                     joggle = TRUE) {
  stopifnot(
    is.list(tensor),
    !is.null(tensor$features_sig),
    !is.null(tensor$features_int),
    is.list(embedding),
    !is.null(embedding$coords_sig),
    !is.null(embedding$coords_int)
  )

  X_sig_3d <- as.matrix(embedding$coords_sig)
  X_int_3d <- as.matrix(embedding$coords_int)

  if (!is.numeric(X_sig_3d)) storage.mode(X_sig_3d) <- "double"
  if (!is.numeric(X_int_3d)) storage.mode(X_int_3d) <- "double"

  if (index < 1L || index > nrow(X_sig_3d) || index > nrow(X_int_3d)) {
    stop("Índice 'index' fora do intervalo do número de circuitries.")
  }

  # barycenters in PCA space (3D)
  bary_sig_3d <- as.numeric(X_sig_3d[index, 1:3])
  bary_int_3d <- as.numeric(X_int_3d[index, 1:3])

  # scale per circuitry from variability in original feature space
  feat_sig <- tensor$features_sig[index, , drop = TRUE]
  feat_int <- tensor$features_int[index, , drop = TRUE]

  sd_sig <- suppressWarnings(stats::sd(as.numeric(feat_sig), na.rm = TRUE))
  sd_int <- suppressWarnings(stats::sd(as.numeric(feat_int), na.rm = TRUE))

  if (!is.finite(sd_sig) || sd_sig <= 0) sd_sig <- 1
  if (!is.finite(sd_int) || sd_int <= 0) sd_int <- 1

  s_sig <- as.numeric(scale_factor) * sd_sig
  s_int <- as.numeric(scale_factor) * sd_int

  make_octahedron <- function(center, s) {
    center <- as.numeric(center)
    s <- as.numeric(s)

    matrix(
      c(
        center[1] + s, center[2],     center[3],
        center[1] - s, center[2],     center[3],
        center[1],     center[2] + s, center[3],
        center[1],     center[2] - s, center[3],
        center[1],     center[2],     center[3] + s,
        center[1],     center[2],     center[3] - s
      ),
      ncol = 3,
      byrow = TRUE,
      dimnames = list(NULL, c("x", "y", "z"))
    )
  }

  pts_sig_3d <- make_octahedron(bary_sig_3d, s_sig)
  pts_int_3d <- make_octahedron(bary_int_3d, s_int)

  # robust hulls
  hull_sig <- .safe_convhull3d(pts_sig_3d, joggle = joggle)
  hull_int <- .safe_convhull3d(pts_int_3d, joggle = joggle)

  list(
    pts_sig_3d  = pts_sig_3d,
    pts_int_3d  = pts_int_3d,
    bary_sig_3d = bary_sig_3d,
    bary_int_3d = bary_int_3d,
    hull_sig    = hull_sig,
    hull_int    = hull_int
  )
}
