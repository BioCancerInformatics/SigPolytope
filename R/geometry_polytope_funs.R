###############################################################################
# geometry_polytope_funs.R
# Funções centrais para geometria latente (tensor 18D → PCA → polytopes 3D)
###############################################################################


###############################################################################
# 0. Nomes fixos das 18 dimensões latentes
###############################################################################

latent_dim_names <- c(
  "rho",              #  1
  "rho_strength",     #  2
  "TN_dir",           #  3
  "TN_strength",      #  4
  "OS_dir",           #  5
  "OS_strength",      #  6
  "OS_lr_chisq",      #  7
  "DSS_dir",          #  8
  "DSS_strength",     #  9
  "DSS_lr_chisq",     # 10
  "DFI_dir",          # 11
  "DFI_strength",     # 12
  "DFI_lr_chisq",     # 13
  "PFI_dir",          # 14
  "PFI_strength",     # 15
  "PFI_lr_chisq",     # 16
  "TME_score",        # 17
  "Immune_dir"        # 18
)

###############################################################################
# 1. Funções auxiliares
###############################################################################

coerce_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(x))
}

neg_log10_p <- function(p, eps = 1e-12) {
  # −log10(p) robusto para NA, zeros, "No data", etc.
  p_num <- coerce_numeric_safe(p)
  out   <- numeric(length(p_num))
  ok    <- !is.na(p_num) & p_num > 0
  out[ok]  <- -log10(p_num[ok] + eps)
  out[!ok] <- 0
  out
}

map_tn_direction <- function(x) {
  # Tumor_vs_normal_*:
  #   "Overexpression"   → +1
  #   "Underexpression"  → −1
  #   "Unchanged"        →  0
  #   "No data" / NA     →  NA (depois viram 0)
  x_chr <- as.character(x)
  x_low <- tolower(trimws(x_chr))
  
  out <- rep(NA_real_, length(x_low))
  out[x_low == "overexpression"]  <-  1
  out[x_low == "underexpression"] <- -1
  out[x_low == "unchanged"]       <-  0
  out
}

map_surv_direction <- function(x) {
  # Cox_*_type_*:
  #   "Risk"       → +1
  #   "Protective" → −1
  #   "NS"/NA      →  0
  x_chr <- as.character(x)
  x_low <- tolower(trimws(x_chr))
  
  out <- rep(0, length(x_low))
  out[grepl("risk",    x_low)] <-  1
  out[grepl("protect", x_low)] <- -1
  out
}

map_immune_direction <- function(x) {
  # Immune_classification_*:
  #   "Immune-hot"/"inflamed"     → +1
  #   "Immune-cold"/"excluded"    → −1
  #   "intermediate/mixed/…"      →  0
  x_chr <- as.character(x)
  x_low <- tolower(trimws(x_chr))
  
  out <- rep(0, length(x_low))
  out[grepl("hot",     x_low)] <-  1
  out[grepl("inflam",  x_low)] <-  1
  out[grepl("cold",    x_low)] <- -1
  out[grepl("exclude", x_low)] <- -1
  out
}

###############################################################################
# 2. Construção do tensor latente 18D a partir de 'circuitries'
###############################################################################

build_geometric_tensor_from_circuitries <- function(circuitries,
                                                    eps = 1e-12,
                                                    verbose = TRUE) {
  stopifnot(is.data.frame(circuitries))
  
  # Colunas obrigatórias (sig/int) — seguem o script original
  required_cols <- c(
    # correlação
    "Correlation_rho_sig", "Correlation_p.adj_sig",
    "Correlation_rho_int", "Correlation_p.adj_int",
    # tumor vs normal (Wilcoxon)
    "Tumor_vs_normal_sig", "Tumor_vs_normal_p.adj_sig",
    "Tumor_vs_normal_int", "Tumor_vs_normal_p.adj_int",
    # OS Cox + log-rank
    "Cox_OS_type_sig",  "Cox_OS_p.value_sig",  "OS_log_rank_chisq_sig",
    "Cox_OS_type_int",  "Cox_OS_p.value_int",  "OS_log_rank_chisq_int",
    # DSS Cox + log-rank
    "Cox_DSS_type_sig", "Cox_DSS_p.value_sig", "DSS_log_rank_chisq_sig",
    "Cox_DSS_type_int", "Cox_DSS_p.value_int", "DSS_log_rank_chisq_int",
    # DFI Cox + log-rank
    "Cox_DFI_type_sig", "Cox_DFI_p.value_sig", "DFI_log_rank_chisq_sig",
    "Cox_DFI_type_int", "Cox_DFI_p.value_int", "DFI_log_rank_chisq_int",
    # PFI Cox + log-rank
    "Cox_PFI_type_sig", "Cox_PFI_p.value_sig", "PFI_log_rank_chisq_sig",
    "Cox_PFI_type_int", "Cox_PFI_p.value_int", "PFI_log_rank_chisq_int",
    # microambiente e imune
    "Microenvironment_score_sig", "Microenvironment_score_int",
    "Immune_classification_sig",  "Immune_classification_int"
  )
  
  missing_cols <- setdiff(required_cols, colnames(circuitries))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in 'circuitries': ",
         paste(missing_cols, collapse = ", "))
  }
  
  N <- nrow(circuitries)
  if (verbose) {
    message("Building 18D geometric tensor from ", N, " circuitries...")
  }
  
  # Helper para um lado (sig OU int)
  build_side_features <- function(d, side = c("sig", "int"), eps = 1e-12) {
    side <- match.arg(side)
    suf  <- if (side == "sig") "_sig" else "_int"
    
    # 1–2) correlação
    rho          <- coerce_numeric_safe(d[[paste0("Correlation_rho", suf)]])
    rho_strength <- neg_log10_p(d[[paste0("Correlation_p.adj", suf)]], eps)
    
    # 3–4) tumor vs normal
    tn_state   <- d[[paste0("Tumor_vs_normal", suf)]]
    tn_dir_raw <- map_tn_direction(tn_state)
    tn_dir     <- tn_dir_raw
    tn_dir[is.na(tn_dir)] <- 0
    
    tn_p        <- d[[paste0("Tumor_vs_normal_p.adj", suf)]]
    tn_p_num    <- coerce_numeric_safe(tn_p)
    tn_strength <- neg_log10_p(tn_p_num, eps)
    
    # 5–16) blocos OS, DSS, DFI, PFI
    surv_names <- c("OS", "DSS", "DFI", "PFI")
    
    surv_dir_list      <- list()
    surv_strength_list <- list()
    surv_chisq_list    <- list()
    
    for (sname in surv_names) {
      type_col <- paste0("Cox_", sname, "_type",   suf)
      p_col    <- paste0("Cox_", sname, "_p.value", suf)
      lr_col   <- paste0(sname, "_log_rank_chisq", suf)
      
      dir  <- map_surv_direction(d[[type_col]])
      dir[is.na(dir)] <- 0
      
      strg <- neg_log10_p(d[[p_col]], eps)
      chis <- coerce_numeric_safe(d[[lr_col]])
      chis[is.na(chis)] <- 0
      
      # zerar força e χ² quando NS
      strg[dir == 0 | is.na(dir)] <- 0
      chis[dir == 0 | is.na(dir)] <- 0
      
      surv_dir_list[[sname]]      <- dir
      surv_strength_list[[sname]] <- strg
      surv_chisq_list[[sname]]    <- chis
    }
    
    # 17) microambiente (score contínuo)
    tme_score <- coerce_numeric_safe(d[[paste0("Microenvironment_score", suf)]])
    tme_score[is.na(tme_score)] <- 0
    
    # 18) direção imune
    immune_dir <- map_immune_direction(d[[paste0("Immune_classification", suf)]])
    immune_dir[is.na(immune_dir)] <- 0
    
    feats <- cbind(
      rho              = rho,
      rho_strength     = rho_strength,
      TN_dir           = tn_dir,
      TN_strength      = tn_strength,
      OS_dir           = surv_dir_list[["OS"]],
      OS_strength      = surv_strength_list[["OS"]],
      OS_lr_chisq      = surv_chisq_list[["OS"]],
      DSS_dir          = surv_dir_list[["DSS"]],
      DSS_strength     = surv_strength_list[["DSS"]],
      DSS_lr_chisq     = surv_chisq_list[["DSS"]],
      DFI_dir          = surv_dir_list[["DFI"]],
      DFI_strength     = surv_strength_list[["DFI"]],
      DFI_lr_chisq     = surv_chisq_list[["DFI"]],
      PFI_dir          = surv_dir_list[["PFI"]],
      PFI_strength     = surv_strength_list[["PFI"]],
      PFI_lr_chisq     = surv_chisq_list[["PFI"]],
      TME_score        = tme_score,
      Immune_dir       = immune_dir
    )
    
    storage.mode(feats) <- "double"
    feats
  }
  
  features_sig <- build_side_features(circuitries, side = "sig", eps = eps)
  features_int <- build_side_features(circuitries, side = "int", eps = eps)
  
  # meta: só rótulos/anotações
  meta_cols <- intersect(
    c(
      "Circuitries_id", "CTAB",
      "Nomenclature_sig", "Nomenclature_int",
      "Signatures", "Interaction",
      "Molecular_class_sig", "Molecular_class_int",
      "Metabolism", "Pathways", "Metabolic_cell_death",
      "Phenotypic_concordance",
      "Immune_concordance",
      "Final_concordance_summary",
      "Microenvironment_classification_sig",
      "Microenvironment_classification_int",
      "Immune_classification_sig",
      "Immune_classification_int",
      "Survival_concordance_OS",
      "Survival_concordance_DSS",
      "Survival_concordance_DFI",
      "Survival_concordance_PFI",
      "Survival_concordance_aggregated",
      "Cox_concordance_OS",
      "Cox_concordance_DSS",
      "Cox_concordance_DFI",
      "Cox_concordance_PFI",
      "Cox_concordance_aggregated",
      "Signature_type", "Interaction_type",
      "Category", "Signature_count", "Interaction_count",
      # 🔹 colunas GEOMÉTRICAS já existentes no seu df:
      "barycenter_distance",
      "distance_implication",
      "sig_hull_vol",
      "int_hull_vol",
      "vol_ratio",
      "vol_implication"
    ),
    colnames(circuitries)
  )
  
  
  meta <- as_tibble(circuitries[, meta_cols, drop = FALSE])
  
  if ("Circuitries_id" %in% colnames(circuitries)) {
    rownames(features_sig) <- circuitries$Circuitries_id
    rownames(features_int) <- circuitries$Circuitries_id
    rownames(meta)         <- circuitries$Circuitries_id
  }
  
  if (verbose) {
    message(
      "Finished building latent feature matrices: ",
      nrow(features_sig), " circuitries × ", ncol(features_sig),
      " dimensions per side (18D)."
    )
  }
  
  list(
    features_sig = features_sig,
    features_int = features_int,
    meta         = meta
  )
}

###############################################################################
# 3. PCA embedding_obj do tensor (18D → 3D)
###############################################################################

compute_pca_embedding_from_tensor <- function(tensor, n_components = 3) {
  
  if (is.null(tensor$features_sig) || is.null(tensor$features_int)) {
    stop("Tensor must contain features_sig and features_int.")
  }
  
  X_sig <- tensor$features_sig
  X_int <- tensor$features_int
  
  if (!is.matrix(X_sig) || !is.matrix(X_int)) {
    stop("tensor$features_sig and tensor$features_int must be matrices.")
  }
  if (ncol(X_sig) != ncol(X_int)) {
    stop("Sig and Int feature matrices must have the same number of columns.")
  }
  
  message("Computing PCA on latent tensor (", 2 * nrow(X_sig),
          " total points × ", ncol(X_sig), " dims)...")
  
  X_all <- rbind(X_sig, X_int)   # 2N × 18
  
  means <- colMeans(X_all, na.rm = TRUE)
  sds   <- apply(X_all, 2, sd, na.rm = TRUE)
  
  zero_var_idx <- which(!is.finite(sds) | sds == 0)
  if (length(zero_var_idx) > 0L) {
    warning("The following latent dimensions have zero variance across all points: ",
            paste(colnames(X_all)[zero_var_idx], collapse = ", "),
            "\nThey will be kept but contribute no variance.")
    sds[zero_var_idx] <- 1
  }
  
  X_scaled <- sweep(X_all, 2, means, "-")
  X_scaled <- sweep(X_scaled, 2, sds,   "/")
  
  pca_obj <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
  pca_obj$center <- means
  pca_obj$scale  <- sds
  
  coords_all <- pca_obj$x[, seq_len(n_components), drop = FALSE]
  
  N <- nrow(X_sig)
  coords_sig <- coords_all[1:N, , drop = FALSE]
  coords_int <- coords_all[(N + 1):(2 * N), , drop = FALSE]
  
  colnames(coords_sig) <- paste0("PC", seq_len(n_components))
  colnames(coords_int) <- paste0("PC", seq_len(n_components))
  
  var_exp <- summary(pca_obj)$importance[2, 1:n_components]
  message(
    "PCA completed. Variance explained by PC1–PC", n_components,
    " = ", round(sum(var_exp) * 100, 2), "%."
  )
  
  list(
    pca        = pca_obj,
    coords_sig = coords_sig,
    coords_int = coords_int
  )
}

###############################################################################
# 4. Projetar vetores latentes arbitrários no espaço PCA
###############################################################################

project_to_pca <- function(pca_obj, X, n_components = 3) {
  X <- as.matrix(X)
  d <- ncol(X)
  
  if (is.null(pca_obj$rotation)) {
    stop("PCA object has no 'rotation'.")
  }
  if (ncol(pca_obj$rotation) != d) {
    stop("Dimension mismatch: X has ", d,
         " columns, PCA rotation has ", ncol(pca_obj$rotation), ".")
  }
  
  X_scaled <- scale(
    X,
    center = pca_obj$center,
    scale  = pca_obj$scale
  )
  scores <- X_scaled %*% pca_obj$rotation
  scores[, seq_len(min(n_components, ncol(scores))), drop = FALSE]
}

###############################################################################
# 5. Construir vértices ± em torno do barycenter (18D)
###############################################################################

build_vertices_for_signature <- function(latent_vec,
                                         side = c("sig", "int"),
                                         pca_obj,
                                         n_components = 3,
                                         eps = 1e-6) {
  side <- match.arg(side)
  v <- as.numeric(latent_vec)
  P <- length(v)   # 18
  
  step <- pmax(abs(v), eps)
  
  verts_latent <- matrix(0, nrow = 2 * P, ncol = P)
  dim_idx      <- integer(2 * P)
  sign_flag    <- character(2 * P)
  
  row_idx <- 1L
  for (j in seq_len(P)) {
    # + step
    vp <- v
    vp[j] <- v[j] + step[j]
    verts_latent[row_idx, ] <- vp
    dim_idx[row_idx]        <- j
    sign_flag[row_idx]      <- "pos"
    row_idx <- row_idx + 1L
    
    # − step
    vn <- v
    vn[j] <- v[j] - step[j]
    verts_latent[row_idx, ] <- vn
    dim_idx[row_idx]        <- j
    sign_flag[row_idx]      <- "neg"
    row_idx <- row_idx + 1L
  }
  
  verts_pca <- project_to_pca(pca_obj, verts_latent, n_components = n_components)
  colnames(verts_pca) <- paste0("PC", seq_len(ncol(verts_pca)))
  
  tibble(
    side      = side,
    dim_index = dim_idx,
    sign      = sign_flag,
    PC1       = verts_pca[, 1],
    PC2       = verts_pca[, 2],
    PC3       = verts_pca[, 3]
  )
}

###############################################################################
# 6. Construir UM polytopo dual (sig + int) para um circuitry
###############################################################################

build_circuitry_polytope <- function(tensor,
                                     embedding_obj,
                                     index = NULL,
                                     circuitry_id = NULL) {
  
  if (is.null(index) && is.null(circuitry_id)) {
    stop("Provide either 'index' or 'circuitry_id'.")
  }
  
  meta <- tensor$meta
  
  # Resolver por Circuitries_id se fornecido
  if (!is.null(circuitry_id)) {
    if (!"Circuitries_id" %in% colnames(meta)) {
      stop("tensor$meta has no 'Circuitries_id' column.")
    }
    idx <- which(meta$Circuitries_id == circuitry_id)
    if (length(idx) == 0L) {
      stop("No circuitry_id '", circuitry_id, "' found in tensor$meta.")
    }
    if (length(idx) > 1L) {
      warning("Multiple rows matched circuitry_id '", circuitry_id,
              "'. Using the first match.")
      idx <- idx[1]
    }
    index <- idx
  }
  
  if (index < 1L || index > nrow(meta)) {
    stop("Index out of range: ", index)
  }
  
  v_sig <- tensor$features_sig[index, ]
  v_int <- tensor$features_int[index, ]
  
  pca_obj <- embedding_obj$pca
  if (is.null(pca_obj)) {
    stop("embedding_obj$pca is NULL.")
  }
  
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
  
  # Rotular dimensão
  verts_sig$dim_name <- latent_dim_names[verts_sig$dim_index]
  verts_int$dim_name <- latent_dim_names[verts_int$dim_index]
  
  vertices_all <- dplyr::bind_rows(verts_sig, verts_int)
  
  # Barycenters projetados
  bary_sig_3d <- project_to_pca(pca_obj, matrix(v_sig, nrow = 1), n_components = 3)[1, ]
  bary_int_3d <- project_to_pca(pca_obj, matrix(v_int, nrow = 1), n_components = 3)[1, ]
  names(bary_sig_3d) <- c("PC1", "PC2", "PC3")
  names(bary_int_3d) <- c("PC1", "PC2", "PC3")
  
  # Hull sig
  hull_sig <- list(faces = NULL, volume = NA_real_, area = NA_real_)
  coords_sig <- as.matrix(verts_sig[, c("PC1", "PC2", "PC3")])
  if (nrow(coords_sig) >= 4L) {
    Hs <- try(geometry::convhulln(coords_sig, options = "FA"), silent = TRUE)
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
  
  # Hull int
  hull_int <- list(faces = NULL, volume = NA_real_, area = NA_real_)
  coords_int <- as.matrix(verts_int[, c("PC1", "PC2", "PC3")])
  if (nrow(coords_int) >= 4L) {
    Hi <- try(geometry::convhulln(coords_int, options = "FA"), silent = TRUE)
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
  
  circuitry_id_out <- if ("Circuitries_id" %in% colnames(meta)) {
    meta$Circuitries_id[index]
  } else {
    as.character(index)
  }
  
  list(
    circuitry_id = circuitry_id_out,
    vertices_sig = verts_sig,
    vertices_int = verts_int,
    vertices_all = vertices_all,
    bary_sig_3d  = bary_sig_3d,
    bary_int_3d  = bary_int_3d,
    hull_sig     = hull_sig,
    hull_int     = hull_int
  )
}

###############################################################################
# 7. Plotly 3D polytope: hulls + vértices + barycenters
###############################################################################

plot_circuitry_polytope <- function(poly_obj,
                                    show_hull = TRUE,
                                    show_vertices = TRUE,
                                    show_barycenters = TRUE,
                                    title = NULL,
                                    outfile_html = NULL) {
  
  verts_all <- poly_obj$vertices_all
  verts_sig <- poly_obj$vertices_sig
  verts_int <- poly_obj$vertices_int
  hull_sig  <- poly_obj$hull_sig
  hull_int  <- poly_obj$hull_int
  bary_sig  <- poly_obj$bary_sig_3d
  bary_int  <- poly_obj$bary_int_3d
  
  if (is.null(title)) {
    title <- paste0("Circuitry polytope (18D → 3D): ", poly_obj$circuitry_id)
  }
  
  p <- plotly::plot_ly()
  
  # 1) Vértices
  if (show_vertices && nrow(verts_all) > 0L) {
    p <- p %>%
      plotly::add_markers(
        data   = verts_all,
        x      = ~PC1, y = ~PC2, z = ~PC3,
        color  = ~side,          # sig vs int
        symbol = ~sign,          # pos vs neg
        symbols = c("circle", "x"),
        marker = list(size = 4),
        text   = ~paste(
          "Side:", side,
          "<br>Dimension:", dim_name,
          "<br>Sign:", sign
        ),
        hoverinfo = "text",
        name  = "Vertices"
      )
  }
  
  # 2) Hull sig
  if (show_hull && !is.null(hull_sig$faces)) {
    simplices <- hull_sig$faces
    if (is.matrix(simplices) && ncol(simplices) == 3) {
      tri   <- as.vector(t(simplices)) - 1L
      n_tri <- nrow(simplices)
      i_idx <- tri[seq(1, 3 * n_tri, by = 3)]
      j_idx <- tri[seq(2, 3 * n_tri, by = 3)]
      k_idx <- tri[seq(3, 3 * n_tri, by = 3)]
      
      coords <- as.matrix(verts_sig[, c("PC1", "PC2", "PC3")])
      
      p <- p %>%
        plotly::add_trace(
          type   = "mesh3d",
          x      = coords[, 1],
          y      = coords[, 2],
          z      = coords[, 3],
          i      = i_idx,
          j      = j_idx,
          k      = k_idx,
          opacity = 0.35,
          name   = "Hull_sig"
        )
    }
  }
  
  # 3) Hull int
  if (show_hull && !is.null(hull_int$faces)) {
    simplices <- hull_int$faces
    if (is.matrix(simplices) && ncol(simplices) == 3) {
      tri   <- as.vector(t(simplices)) - 1L
      n_tri <- nrow(simplices)
      i_idx <- tri[seq(1, 3 * n_tri, by = 3)]
      j_idx <- tri[seq(2, 3 * n_tri, by = 3)]
      k_idx <- tri[seq(3, 3 * n_tri, by = 3)]
      
      coords <- as.matrix(verts_int[, c("PC1", "PC2", "PC3")])
      
      p <- p %>%
        plotly::add_trace(
          type   = "mesh3d",
          x      = coords[, 1],
          y      = coords[, 2],
          z      = coords[, 3],
          i      = i_idx,
          j      = j_idx,
          k      = k_idx,
          opacity = 0.35,
          name   = "Hull_int"
        )
    }
  }
  
  # 4) Barycenters
  if (show_barycenters) {
    p <- p %>%
      plotly::add_markers(
        x = bary_sig["PC1"],
        y = bary_sig["PC2"],
        z = bary_sig["PC3"],
        marker   = list(size = 8, symbol = "diamond"),
        name     = "Barycenter_sig",
        text     = "Signature barycenter (18D → 3D)",
        hoverinfo = "text"
      ) %>%
      plotly::add_markers(
        x = bary_int["PC1"],
        y = bary_int["PC2"],
        z = bary_int["PC3"],
        marker   = list(size = 8, symbol = "diamond-open"),
        name     = "Barycenter_int",
        text     = "Interaction barycenter (18D → 3D)",
        hoverinfo = "text"
      )
  }
  
  # 5) Legend com as 18 dimensões
  dim_legend_text <- paste0(
    "<b>Latent dimensions (18D):</b><br>",
    paste0(seq_along(latent_dim_names), " = ", latent_dim_names, collapse = "<br>")
  )
  
  p <- p %>%
    plotly::layout(
      title = list(text = title),
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3")
      ),
      legend = list(
        x = 0.98, y = 1.0,
        xanchor = "right",
        yanchor = "top"
      ),
      annotations = list(
        list(
          x = 0.02, y = 0.98,
          xref = "paper", yref = "paper",
          showarrow = FALSE,
          align = "left",
          text = dim_legend_text,
          bordercolor = "black",
          borderwidth = 0.5,
          bgcolor = "rgba(255,255,255,0.7)"
        )
      )
    )
  
  if (!is.null(outfile_html)) {
    safe_file <- gsub("[^A-Za-z0-9_\\.]+", "_", outfile_html)
    htmlwidgets::saveWidget(plotly::as_widget(p), file = safe_file, selfcontained = TRUE)
    message("Polytope HTML saved to: ", normalizePath(safe_file))
  }
  
  p
}

###############################################################################
# 8. Distância de barycenters = discordância multidimensional
###############################################################################

barycenter_distance <- function(i, tensor, embedding_obj) {
  poly_i <- build_circuitry_polytope(tensor, embedding_obj, index = i)
  sqrt(sum((poly_i$bary_sig_3d - poly_i$bary_int_3d)^2))
}



###############################################################################
# 9. add_geometry_metadata — versão leve (usar valores já existentes no df)
###############################################################################

add_geometry_metadata <- function(tensor, embedding_obj, ...) {
  # 👉 Como as colunas geométricas já existem no seu data.frame original
  # (barycenter_distance, distance_implication, sig_hull_vol,
  #  int_hull_vol, vol_ratio, vol_implication),
  # não precisamos recalcular nada aqui.
  #
  # Esta função vira apenas um "pass-through" para manter compatibilidade.
  tensor
}
