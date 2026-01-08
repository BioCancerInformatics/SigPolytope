# =========================================================
# plot_circuitry_3d.R  (1 engine + 3 wrappers + fixed palette)
# =========================================================

# -----------------------------
# Default palette (matches your screenshot style)
# -----------------------------
# Signature = teal hull + dark blue points/bary
# Regulator = orange hull + dark orange points/bary
#' @keywords internal
GEOCIRCUITRY_PAL <- list(
  sig_faces = "#56B4E9",
  sig_pts   = "#0072B2",
  sig_bary  = "#0072B2",
  int_faces = "#E69F00",
  int_pts   = "#D55E00",
  int_bary  = "#D55E00",
  dist_line = "#000000"
)

# -----------------------------
# Internal helper: null-coalesce
# -----------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# Internal helper: extract faces robustly (new + legacy schemas)
# -----------------------------
#' @keywords internal
.extract_hull_faces <- function(hull) {
  if (is.null(hull)) return(NULL)

  # legacy: hull$hull$hull
  if (is.list(hull) && !is.null(hull$hull) && is.list(hull$hull) && !is.null(hull$hull$hull)) {
    return(hull$hull$hull)
  }

  # new: hull$faces
  if (is.list(hull) && !is.null(hull$faces) && is.matrix(hull$faces)) {
    return(hull$faces)
  }

  # direct matrix
  if (is.matrix(hull)) return(hull)

  NULL
}

# -----------------------------
# Internal helper: validate palette keys
# -----------------------------
#' @keywords internal
.validate_palette <- function(palette) {
  req <- c("sig_faces","sig_pts","sig_bary","int_faces","int_pts","int_bary","dist_line")
  miss <- setdiff(req, names(palette))
  if (length(miss) > 0) {
    stop("Palette inválida. Faltam chaves: ", paste(miss, collapse = ", "))
  }
  palette
}

# -----------------------------
# Internal engine (do NOT export)
# -----------------------------
#' @keywords internal
.plot_poly_3d <- function(points_3d,
                          bary_3d = NULL,
                          hull = NULL,
                          title = NULL,
                          alpha_faces = 0.25,
                          size_vertices = 4,
                          size_bary = 7,
                          color_vertices = "#0072B2",
                          color_faces    = "#56B4E9",
                          color_bary     = "#0072B2",
                          symbol_vertices = "circle",
                          symbol_bary = "circle",
                          showlegend_vertices = TRUE,
                          showlegend_bary = TRUE,
                          name_vertices = "Vertices",
                          name_bary = "Barycenter",
                          legendgroup = "poly",
                          dims_text = NULL) {

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("O pacote 'plotly' é necessário para plotar em 3D.")
  }

  pts <- as.data.frame(points_3d)
  colnames(pts) <- c("x", "y", "z")

  p <- plotly::plot_ly()

  # hull (para ficar "atrás" dos pontos visualmente)
  faces <- .extract_hull_faces(hull)
  if (!is.null(faces)) {

    # força uma cor por face (mais robusto que 'color=')
    face_cols <- rep(color_faces, nrow(faces))

    p <- p |>
      plotly::add_trace(
        type = "mesh3d",
        x = pts$x, y = pts$y, z = pts$z,
        i = faces[, 1] - 1L,
        j = faces[, 2] - 1L,
        k = faces[, 3] - 1L,
        name = paste0(name_vertices, " hull"),
        legendgroup = legendgroup,
        showlegend  = FALSE,
        opacity     = alpha_faces,
        facecolor   = face_cols,
        flatshading = TRUE
      )
  }

  # vertices
  p <- p |>
    plotly::add_markers(
      data = pts,
      x = ~x, y = ~y, z = ~z,
      name = name_vertices,
      legendgroup = legendgroup,
      showlegend = showlegend_vertices,
      marker = list(
        color  = color_vertices,
        size   = size_vertices,
        symbol = symbol_vertices
      )
    )

  # barycenter
  if (!is.null(bary_3d)) {
    b <- setNames(as.numeric(bary_3d), c("x", "y", "z"))
    p <- p |>
      plotly::add_markers(
        x = b["x"], y = b["y"], z = b["z"],
        name = name_bary,
        legendgroup = legendgroup,
        showlegend = showlegend_bary,
        marker = list(
          color  = color_bary,
          size   = size_bary,
          symbol = symbol_bary
        )
      )
  }

  ann <- NULL
  if (!is.null(dims_text)) {
    ann <- list(
      list(
        text = dims_text,
        x = 0,
        y = 1,
        xref = "paper",
        yref = "paper",
        xanchor = "left",
        yanchor = "top",
        align = "left",
        showarrow = FALSE,
        font = list(size = 11),
        bgcolor = "rgba(255,255,255,0.85)",
        bordercolor = "black",
        borderwidth = 1
      )
    )
  }


  # layout
  p <- p |>
    plotly::layout(
      title = title,
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3")
      ),
      legend = list(orientation = "h", x = 0, y = 1.05),
      annotations = ann
    )

  p
}

# -----------------------------
# Helper: accept either (sig_geom, reg_geom, index) or a ready poly object
# -----------------------------
#' @keywords internal
.resolve_poly_input <- function(poly = NULL, sig_geom = NULL, reg_geom = NULL, index = 1) {
  if (!is.null(poly)) return(poly)

  if (is.null(sig_geom) || is.null(reg_geom)) {
    stop("Forneça `poly` OU (`sig_geom` e `reg_geom`).")
  }

  aligned <- .align_sig_int_features(sig_geom, reg_geom, strict = FALSE)

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

  embedding <- compute_pca_embedding_from_tensor(tensor, n_components = 3)
  build_circuitry_polytope(tensor, embedding, index = index)
}


# =========================================================
# Public wrappers (palette is fixed internally by default)
# =========================================================

#' 3D polytope da ASSINATURA (signature) para um circuitry
#'
#' @inheritParams .plot_poly_3d
#' @return plotly
#' @export
plot_signature_polytope_3d <- function(poly = NULL,
                                       sig_geom = NULL,
                                       reg_geom = NULL,
                                       index = 1,
                                       alpha_faces   = 0.25,
                                       size_vertices = 4,
                                       size_bary     = 7) {

  palette <- .validate_palette(GEOCIRCUITRY_PAL)
  poly <- .resolve_poly_input(poly, sig_geom, reg_geom, index)

  # -----------------------------
  # Annotation: latent dimensions used in the geometry
  # -----------------------------
  dims_text <- NULL
  if (!is.null(poly$latent_dims) && !is.null(poly$latent_dims$signature)) {
    dims_text <- .format_latent_dims_text(
      poly$latent_dims$signature,
      title = "Latent dimensions (signature)"
    )
  } else if (!is.null(poly$latent_dims) && !is.null(poly$latent_dims$sig)) {
    dims_text <- .format_latent_dims_text(
      poly$latent_dims$sig,
      title = "Latent dimensions (signature)"
    )
  } else if (!is.null(poly$dim_names_sig)) {
    dims_text <- .format_latent_dims_text(
      poly$dim_names_sig,
      title = "Latent dimensions (signature)"
    )
  }

  ann <- NULL
  if (!is.null(dims_text)) {
    ann <- list(
      list(
        text = dims_text,
        x = 0,
        y = 1,
        xref = "paper",
        yref = "paper",
        xanchor = "left",
        yanchor = "top",
        align = "left",
        showarrow = FALSE,
        font = list(size = 11),
        bgcolor = "rgba(255,255,255,0.85)",
        bordercolor = "black",
        borderwidth = 1
      )
    )
  }




  ###


  .plot_poly_3d(
    points_3d = poly$pts_sig_3d,
    bary_3d   = poly$bary_sig_3d,
    hull      = poly$hull_sig,
    title     = "Signature polytope",
    alpha_faces   = alpha_faces,
    size_vertices = size_vertices,
    size_bary     = size_bary,
    color_vertices = palette$sig_pts,
    color_faces    = palette$sig_faces,
    color_bary     = palette$sig_bary,
    symbol_vertices = "circle",
    symbol_bary     = "square",
    name_vertices = "Signature vertices",
    name_bary     = "Signature barycenter",
    legendgroup   = "sig",
    dims_text     = dims_text
  )
}

#' 3D polytope do REGULADOR para um circuitry
#'
#' @inheritParams .plot_poly_3d
#' @return plotly
#' @export
plot_regulator_polytope_3d <- function(poly = NULL,
                                       sig_geom = NULL,
                                       reg_geom = NULL,
                                       index = 1,
                                       alpha_faces   = 0.25,
                                       size_vertices = 4,
                                       size_bary     = 7) {

  palette <- .validate_palette(GEOCIRCUITRY_PAL)
  poly <- .resolve_poly_input(poly, sig_geom, reg_geom, index)

  # -----------------------------
  # Annotation: latent dimensions (regulator)
  # -----------------------------
  dims_text <- NULL

  if (!is.null(poly$latent_dims) && !is.null(poly$latent_dims$regulator)) {
    dims_text <- .format_latent_dims_text(
      poly$latent_dims$regulator,
      title = "Latent dimensions (regulator)"
    )
  } else if (!is.null(poly$latent_dims) && !is.null(poly$latent_dims$int)) {
    dims_text <- .format_latent_dims_text(
      poly$latent_dims$int,
      title = "Latent dimensions (regulator)"
    )
  } else if (!is.null(poly$dim_names_int)) {
    dims_text <- .format_latent_dims_text(
      poly$dim_names_int,
      title = "Latent dimensions (regulator)"
    )
  }


  .plot_poly_3d(
    points_3d = poly$pts_int_3d,
    bary_3d   = poly$bary_int_3d,
    hull      = poly$hull_int,
    title     = "Regulator polytope",
    alpha_faces   = alpha_faces,
    size_vertices = size_vertices,
    size_bary     = size_bary,
    color_vertices = palette$int_pts,
    color_faces    = palette$int_faces,
    color_bary     = palette$int_bary,
    symbol_vertices = "diamond",
    symbol_bary     = "square",
    name_vertices = "Regulator vertices",
    name_bary     = "Regulator barycenter",
    legendgroup   = "int",
    dims_text     = dims_text
  )
}

#' 3D circuitry convergence (signature x regulator)
#'
#' Plota assinatura + regulador + linha barycenter→barycenter,
#' e preenche os hulls com as cores padrão do pacote.
#'
#' @param alpha_faces Transparência dos hulls.
#' @param size_vertices Tamanho dos vértices.
#' @param size_bary Tamanho dos baricentros.
#' @return plotly
#' @export
plot_circuitry_convergence_3d <- function(
    poly = NULL,
    sig_geom = NULL,
    reg_geom = NULL,
    index = 1,
    alpha_faces   = 0.25,
    size_vertices = 4,
    size_bary     = 7
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("O pacote 'plotly' é necessário para plot_circuitry_convergence_3d().")
  }

  palette <- .validate_palette(GEOCIRCUITRY_PAL)
  poly <- .resolve_poly_input(poly, sig_geom, reg_geom, index)
  # -----------------------------
  # Annotation: latent dimensions (regulator)
  # -----------------------------
  dims_text <- NULL

  if (!is.null(poly$latent_dims) && !is.null(poly$latent_dims$regulator)) {
    dims_text <- .format_latent_dims_text(
      poly$latent_dims$regulator,
      title = "Latent dimensions (regulator)"
    )
  } else if (!is.null(poly$latent_dims) && !is.null(poly$latent_dims$int)) {
    dims_text <- .format_latent_dims_text(
      poly$latent_dims$int,
      title = "Latent dimensions (regulator)"
    )
  } else if (!is.null(poly$dim_names_int)) {
    dims_text <- .format_latent_dims_text(
      poly$dim_names_int,
      title = "Latent dimensions (regulator)"
    )
  }


  sig_df <- as.data.frame(poly$pts_sig_3d); colnames(sig_df) <- c("x","y","z")
  int_df <- as.data.frame(poly$pts_int_3d); colnames(int_df) <- c("x","y","z")

  bary_s <- setNames(as.numeric(poly$bary_sig_3d), c("x","y","z"))
  bary_i <- setNames(as.numeric(poly$bary_int_3d), c("x","y","z"))
  d_bary <- sqrt(sum((bary_s - bary_i)^2))

  p <- plotly::plot_ly()

  # hulls (faces) first
  faces_sig <- .extract_hull_faces(poly$hull_sig)
  if (!is.null(faces_sig)) {
    p <- p |>
      plotly::add_trace(
        type = "mesh3d",
        x = sig_df$x, y = sig_df$y, z = sig_df$z,
        i = faces_sig[,1] - 1L,
        j = faces_sig[,2] - 1L,
        k = faces_sig[,3] - 1L,
        name = "Signature hull",
        legendgroup = "sig",
        opacity = alpha_faces,
        color   = palette$sig_faces,
        showlegend = FALSE
      )
  }

  faces_int <- .extract_hull_faces(poly$hull_int)
  if (!is.null(faces_int)) {
    p <- p |>
      plotly::add_trace(
        type = "mesh3d",
        x = int_df$x, y = int_df$y, z = int_df$z,
        i = faces_int[,1] - 1L,
        j = faces_int[,2] - 1L,
        k = faces_int[,3] - 1L,
        name = "Regulator hull",
        legendgroup = "int",
        opacity = alpha_faces,
        color   = palette$int_faces,
        showlegend = FALSE
      )
  }

  # vertices + bary (after hull)
  p <- p |>
    plotly::add_markers(
      data = sig_df, x = ~x, y = ~y, z = ~z,
      name = "Signature vertices",
      legendgroup = "sig",
      marker = list(color = palette$sig_pts, size = size_vertices, symbol = "circle")
    ) |>
    plotly::add_markers(
      x = bary_s["x"], y = bary_s["y"], z = bary_s["z"],
      name = "Signature barycenter",
      legendgroup = "sig",
      marker = list(color = palette$sig_bary, size = size_bary, symbol = "square")
    )

  p <- p |>
    plotly::add_markers(
      data = int_df, x = ~x, y = ~y, z = ~z,
      name = "Regulator vertices",
      legendgroup = "int",
      marker = list(color = palette$int_pts, size = size_vertices, symbol = "diamond")
    ) |>
    plotly::add_markers(
      x = bary_i["x"], y = bary_i["y"], z = bary_i["z"],
      name = "Regulator barycenter",
      legendgroup = "int",
      marker = list(color = palette$int_bary, size = size_bary, symbol = "square")
    )

  # bary→bary line
  p <- p |>
    plotly::add_trace(
      x = c(bary_s["x"], bary_i["x"]),
      y = c(bary_s["y"], bary_i["y"]),
      z = c(bary_s["z"], bary_i["z"]),
      type = "scatter3d",
      mode = "lines",
      name = "Barycenter distance",
      legendgroup = "dist",
      line = list(color = palette$dist_line, width = 4, dash = "dash")
    )

  circuitry_id <- poly$circuitry_id %||% paste0("Circuitry_", index)

  p |>
    plotly::layout(
      title = paste0("Circuitry ", circuitry_id,
                     " — barycenter distance = ", sprintf("%.2f", d_bary)),
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3")
      ),
      legend = list(orientation = "h", x = 0, y = 1.05)
    )
}

#' 3D circuitry polytope (signature x regulator) with hulls
#'
#' Wrapper de compatibilidade: equivalente a plot_circuitry_convergence_3d().
#'
#' @inheritParams plot_circuitry_convergence_3d
#' @return plotly
#' @export
plot_circuitry_3d <- function(
    poly = NULL,
    sig_geom = NULL,
    reg_geom = NULL,
    index = 1,
    alpha_faces   = 0.25,
    size_vertices = 4,
    size_bary     = 7
) {
  plot_circuitry_convergence_3d(
    poly = poly,
    sig_geom = sig_geom,
    reg_geom = reg_geom,
    index = index,
    alpha_faces = alpha_faces,
    size_vertices = size_vertices,
    size_bary = size_bary
  )
}
