# Em build_circuitry_polytope.R - RENOMEAR para:
#' Extract basic geometry for a single circuitry
#'
#' Função mais simples que extrai apenas a geometria básica
#' (signature + regulator) para um circuitry.
#'
#' @param tensor Lista com `features_sig`, `features_int` e `meta`
#' @param embedding Lista com `coords_sig`, `coords_int` e `pca`
#' @param index Índice (linha) da circuitry na tabela (opcional se usar `circuitry_id`)
#' @param circuitry_id ID da circuitry (ex.: "LGG-3207 / LGG-2427")
#'
#' @export
extract_circuitry_geometry <- function(tensor, embedding,
                                       index = NULL,
                                       circuitry_id = NULL) {

  meta <- tensor$meta

  # escolher a linha-alvo
  if (!is.null(circuitry_id)) {
    if (!("Circuitries_id" %in% colnames(meta))) {
      stop("Coluna 'Circuitries_id' não encontrada em tensor$meta.")
    }
    idx <- match(circuitry_id, meta$Circuitries_id)
    if (is.na(idx)) {
      stop("circuitry_id não encontrado em tensor$meta$Circuitries_id: ", circuitry_id)
    }
  } else if (!is.null(index)) {
    idx <- as.integer(index)
    if (idx < 1 || idx > nrow(meta)) {
      stop("index fora do intervalo válido.")
    }
  } else {
    stop("Forneça ou `index` ou `circuitry_id`.")
  }

  cid <- meta$Circuitries_id[idx]

  sig_geom <- list(
    side   = "sig",
    id     = cid,
    latent = tensor$features_sig[idx, , drop = FALSE],
    coords = embedding$coords_sig[idx, , drop = FALSE],
    pca    = embedding$pca,
    meta   = meta[idx, , drop = FALSE]
  )

  int_geom <- list(
    side   = "int",
    id     = cid,
    latent = tensor$features_int[idx, , drop = FALSE],
    coords = embedding$coords_int[idx, , drop = FALSE],
    pca    = embedding$pca,
    meta   = meta[idx, , drop = FALSE]
  )

  list(
    circuitry_id = cid,
    sig_geom     = sig_geom,
    int_geom     = int_geom
  )
}



#' Internal: build preview octahedron polytope (NOT paper method)
#'
#' Cria um octaedro 3D simples ao redor dos baricentros (em PCA 3D),
#' com escala baseada na variabilidade das features (por circuitry).
#' Função de *preview rápido* (não é o método do paper).
#'
#' Melhorias (robustez):
#' - usa MAD como fallback do SD (mais robusto)
#' - impõe escala mínima (min_scale) para evitar colapso visual/numérico
#' - se a nuvem ficar degenerada (rank < 3), aplica jitter adaptativo e tenta hull com QJ
#'
#' @param tensor Lista com features_sig, features_int, meta
#' @param embedding Lista com coords_sig, coords_int (PCA 3D)
#' @param index Índice do circuitry (1..N)
#' @param scale_factor Fator de escala
#' @param min_scale Escala mínima (em unidades de PCA) para evitar politopo colapsado
#' @param use_mad Se TRUE, usa MAD como medida principal de variabilidade (fallback para SD)
#' @param jitter_frac Fração do "span" da nuvem usada para jitter quando rank < 3
#'
#' @return Lista com vértices e hulls (preview)
#' @keywords internal
build_preview_octahedron <- function(tensor,
                                     embedding,
                                     index,
                                     scale_factor = 0.5,
                                     min_scale    = 0.05,
                                     use_mad      = TRUE,
                                     jitter_frac  = 1e-6) {
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

  if (nrow(X_sig_3d) < index || nrow(X_int_3d) < index) {
    stop("Índice 'index' fora do intervalo do número de circuitries.")
  }

  # barycenters em PCA 3D
  bary_sig_3d <- as.numeric(X_sig_3d[index, 1:3, drop = TRUE])
  bary_int_3d <- as.numeric(X_int_3d[index, 1:3, drop = TRUE])

  # features originais (18D ou P dims) do circuitry
  feat_sig <- tensor$features_sig[index, , drop = TRUE]
  feat_int <- tensor$features_int[index, , drop = TRUE]

  feat_sig <- suppressWarnings(as.numeric(feat_sig))
  feat_int <- suppressWarnings(as.numeric(feat_int))
  feat_sig[!is.finite(feat_sig)] <- 0
  feat_int[!is.finite(feat_int)] <- 0

  # medida de variabilidade (robusta)
  spread1 <- function(v) {
    v <- as.numeric(v)
    v <- v[is.finite(v)]
    if (length(v) < 2) return(NA_real_)

    if (isTRUE(use_mad)) {
      m <- stats::mad(v, constant = 1, na.rm = TRUE)
      if (is.finite(m) && m > 0) return(m)
    }

    s <- stats::sd(v, na.rm = TRUE)
    if (is.finite(s) && s > 0) return(s)

    # fallback final: amplitude
    r <- diff(range(v, na.rm = TRUE))
    if (is.finite(r) && r > 0) return(r / 2)

    NA_real_
  }

  sp_sig <- spread1(feat_sig)
  sp_int <- spread1(feat_int)

  # escala final (força mínimo)
  s_sig <- max(min_scale, scale_factor * ifelse(is.finite(sp_sig), sp_sig, 1))
  s_int <- max(min_scale, scale_factor * ifelse(is.finite(sp_int), sp_int, 1))

  make_octahedron <- function(center, s) {
    center <- as.numeric(center)
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

  # hull robusto (com rank check + jitter/QJ quando necessário)
  hs <- .safe_convhull3d(pts_sig_3d, joggle = TRUE, jitter_frac = jitter_frac)
  hi <- .safe_convhull3d(pts_int_3d, joggle = TRUE, jitter_frac = jitter_frac)

  list(
    pts_sig_3d  = pts_sig_3d,
    pts_int_3d  = pts_int_3d,
    bary_sig_3d = bary_sig_3d,
    bary_int_3d = bary_int_3d,
    hull_sig    = hs,
    hull_int    = hi,
    scale_sig   = s_sig,
    scale_int   = s_int
  )
}


#' @keywords internal
.safe_convhull3d <- function(coords,
                             joggle = TRUE,
                             jitter_frac = 1e-6,
                             max_tries = 3) {
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("Precisa do pacote 'geometry' para convhulln().")
  }

  coords <- as.matrix(coords)
  coords <- coords[stats::complete.cases(coords), , drop = FALSE]
  coords <- unique(coords)

  if (nrow(coords) < 4) {
    return(list(faces = NULL, volume = 0, area = 0))
  }

  # rank afim (se < 3, não há volume 3D)
  rank_affine <- function(X) {
    Xc <- scale(X, center = TRUE, scale = FALSE)
    s <- svd(Xc)$d
    if (!length(s)) return(0L)
    sum(s > max(s) * 1e-12)
  }

  r <- rank_affine(coords)

  # jitter adaptativo (escala pelo "span" da nuvem)
  add_jitter <- function(X, frac) {
    span <- apply(X, 2, function(v) diff(range(v)))
    span[!is.finite(span) | span == 0] <- 1
    eps <- frac * max(span)
    X + matrix(stats::runif(length(X), -eps, eps), nrow = nrow(X))
  }

  # tenta convex hull com escalonamento de robustez
  try_hull <- function(X, use_qj) {
    opt <- if (use_qj) "FA QJ" else "FA"
    H <- try(geometry::convhulln(X, options = opt), silent = TRUE)
    if (inherits(H, "try-error")) return(NULL)

    faces <- H
    vol   <- attr(H, "vol")
    area  <- attr(H, "area")

    # algumas versões retornam lista
    if (is.list(H)) {
      faces <- if (!is.null(H$hull)) H$hull else if (!is.null(H$facets)) H$facets else NULL
      if (!is.null(H$vol))  vol  <- H$vol
      if (!is.null(H$area)) area <- H$area
    }

    list(faces = faces, volume = as.numeric(vol), area = as.numeric(area))
  }

  # Se rank < 3, tentamos “desdegenerar” com jitter pequeno (somente para robustez numérica/preview)
  Xwork <- coords
  if (r < 3) {
    # 1) tenta com QJ direto (às vezes resolve)
    out <- try_hull(Xwork, use_qj = TRUE)
    if (!is.null(out)) return(out)

    # 2) jitter progressivo
    for (k in seq_len(max_tries)) {
      Xj <- add_jitter(Xwork, frac = jitter_frac * (10^k))
      if (rank_affine(Xj) >= 3) {
        out <- try_hull(Xj, use_qj = TRUE)
        if (!is.null(out)) return(out)
      }
    }

    # ainda degenerado: sem volume
    return(list(faces = NULL, volume = 0, area = 0))
  }

  # rank OK: tenta normal, depois QJ se falhar
  out <- try_hull(Xwork, use_qj = FALSE)
  if (!is.null(out)) return(out)

  if (isTRUE(joggle)) {
    out2 <- try_hull(Xwork, use_qj = TRUE)
    if (!is.null(out2)) return(out2)
  }

  # última tentativa: jitter leve + QJ
  for (k in seq_len(max_tries)) {
    Xj <- add_jitter(Xwork, frac = jitter_frac * (10^k))
    out3 <- try_hull(Xj, use_qj = TRUE)
    if (!is.null(out3)) return(out3)
  }

  list(faces = NULL, volume = 0, area = 0)
}
