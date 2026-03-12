## ===============================
## NixOS / RStudio Server HARDENING
## ===============================
if (nzchar(Sys.which("chromium"))) {
  tmp_profile <- tempfile("chromote-profile-")
  dir.create(tmp_profile)
  
  Sys.setenv(
    CHROMOTE_CHROME = Sys.which("chromium"),
    CHROMOTE_ARGS   = paste(
      "--no-sandbox",
      "--disable-dev-shm-usage",
      "--headless=new",
      "--remote-debugging-port=0",
      sprintf("--user-data-dir=%s", tmp_profile)
    )
  )
}

###############################################################################
# Geometric Multidimensional Representation of Regulatory Circuitries
# -------------------------------------------------------------------
# Notes:
# (1) The barycenter distance is the Euclidean distance between the two 3D barycenters in PCA space, acting as your “multidimensional discordance magnitude.”
# (2) In geometric terms the barycenter-distance behaves exactly like a multidimensional “discordance magnitude.”
# In PCA units:

# Distance < 0.5 → high concordance
# 
# Distance 0.5–1.5 → moderate concordance or mild discordance
# 
# Distance > 1.5 → strong discordance
# 
# Distance > 2.5 → extreme discordance
#  
# 18D latent tensor per side (signature / interaction) + dual polytope hulls
#
# This script builds:
#   1) An 18D latent tensor from the 'circuitries' table, for each side:
#        - sig  (signature)
#        - int  (interaction)
#   2) A global PCA embedding (3D) on all 2N latent vectors.
#   3) For each circuitry (row), TWO polytopes in the same 3D PCA space:
#        - one convex hull for the signature side (sig),
#        - one convex hull for the interaction side (int),
#      each built from 18 barycenter-centred ± vertices.
#   4) A Plotly 3D visualization with:
#        - sig vs int hulls,
#        - sig/int vertices,
#        - sig/int barycenters,
#        - an annotation describing the 18 latent dimensions.
#
# Latent schema per side (18D):
#   1)  rho             : correlation coefficient                 (Correlation_rho_*)
#   2)  rho_strength    : −log10(adjusted p for correlation)      (Correlation_p.adj_*)
#   3)  TN_dir          : tumor vs normal direction               (Tumor_vs_normal_*)
#   4)  TN_strength     : −log10(Wilcoxon p)                      (Tumor_vs_normal_p.adj_*)
#   5)  OS_dir          : Cox OS direction                        (Cox_OS_type_*)
#   6)  OS_strength     : −log10(Cox OS p)                        (Cox_OS_p.value_*)
#   7)  OS_lr_chisq     : log-rank χ² for OS                      (OS_log_rank_chisq_*)
#   8)  DSS_dir         : Cox DSS direction                       (Cox_DSS_type_*)
#   9)  DSS_strength    : −log10(Cox DSS p)                       (Cox_DSS_p.value_*)
#   10) DSS_lr_chisq    : log-rank χ² for DSS                     (DSS_log_rank_chisq_*)
#   11) DFI_dir         : Cox DFI direction                       (Cox_DFI_type_*)
#   12) DFI_strength    : −log10(Cox DFI p)                       (Cox_DFI_p.value_*)
#   13) DFI_lr_chisq    : log-rank χ² for DFI                     (DFI_log_rank_chisq_*)
#   14) PFI_dir         : Cox PFI direction                       (Cox_PFI_type_*)
#   15) PFI_strength    : −log10(Cox PFI p)                       (Cox_PFI_p.value_*)
#   16) PFI_lr_chisq    : log-rank χ² for PFI                     (PFI_log_rank_chisq_*)
#   17) TME_score       : microenvironment score                  (Microenvironment_score_*)
#   18) Immune_dir      : immune class direction                  (Immune_classification_*)
#
# Important survival rule:
#   For each endpoint (OS, DSS, DFI, PFI), if Cox_*_type_* == "NS"
#   (or maps to direction 0), then:
#        dir      = 0
#        strength = 0
#        lr_chisq = 0
#   → NS endpoints do not deform the polytope along those axes.
###############################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("dplyr", quietly = TRUE))        install.packages("dplyr")
  if (!requireNamespace("tibble", quietly = TRUE))       install.packages("tibble")
  if (!requireNamespace("geometry", quietly = TRUE))     install.packages("geometry")
  if (!requireNamespace("plotly", quietly = TRUE))       install.packages("plotly")
  if (!requireNamespace("htmlwidgets", quietly = TRUE))  install.packages("htmlwidgets")
})

library(rio)
library(dplyr)
library(tibble)
library(geometry)
library(plotly)
library(htmlwidgets)
library(parallel)
library(ggplot2)
library(readr)

###############################################################################
# 0. OPTIONAL: import circuitries from TSV if not already in memory
###############################################################################
# If `circuitries` is already loaded, skip this.
if (!requireNamespace("rio", quietly = TRUE)) install.packages("rio")
library(rio)
circuitries <- import("Regulatory_circuitries_02.tsv")

###############################################################################
# 0. Latent dimension names (18D)
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
# 1. Helper functions for encoding and numeric safety
###############################################################################

coerce_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(x))
}

neg_log10_p <- function(p, eps = 1e-12) {
  # Generic −log10(p) with robust handling of "No data", NA, zeros, etc.
  p_num <- coerce_numeric_safe(p)
  out   <- numeric(length(p_num))
  ok    <- !is.na(p_num) & p_num > 0
  out[ok]  <- -log10(p_num[ok] + eps)
  out[!ok] <- 0
  out
}

map_tn_direction <- function(x) {
  # Tumor_vs_normal_* categories:
  #   "Overexpression"   → +1
  #   "Underexpression"  → −1
  #   "Unchanged"        →  0
  #   "No data" / NA     →  NA (later set to 0)
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
  #   "NS" / NA    →  0
  x_chr <- as.character(x)
  x_low <- tolower(trimws(x_chr))
  
  out <- rep(0, length(x_low))
  out[grepl("risk",    x_low)] <-  1
  out[grepl("protect", x_low)] <- -1
  out
}

map_immune_direction <- function(x) {
  # Immune_classification_*:
  #   "Immune-hot" / "inflamed"      → +1
  #   "Immune-cold" / "excluded"     → −1
  #   "intermediate", "mixed", etc.  →  0
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
# 2. Build 18D latent tensor from 'circuitries'
###############################################################################
# Per side (sig / int), 18 coordinates with four survival endpoints (OS, DSS,
# DFI, PFI). Tumor–normal p-values are Wilcoxon-based (not FDR), even if the
# column is labeled "*_p.adj".
###############################################################################

build_geometric_tensor_from_circuitries <- function(circuitries,
                                                    eps = 1e-12,
                                                    verbose = TRUE) {
  stopifnot(is.data.frame(circuitries))
  
  # Required columns for all endpoints
  required_cols <- c(
    # correlation
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
    # microenvironment and immune
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
  
  # Helper to build features for one side
  build_side_features <- function(d, side = c("sig", "int"), eps = 1e-12) {
    side <- match.arg(side)
    suf  <- if (side == "sig") "_sig" else "_int"
    
    # 1–2) Correlation
    rho          <- coerce_numeric_safe(d[[paste0("Correlation_rho", suf)]])
    rho_strength <- neg_log10_p(d[[paste0("Correlation_p.adj", suf)]], eps)
    
    # 3–4) Tumor vs normal
    tn_state   <- d[[paste0("Tumor_vs_normal", suf)]]
    tn_dir_raw <- map_tn_direction(tn_state)
    tn_dir     <- tn_dir_raw
    tn_dir[is.na(tn_dir)] <- 0
    
    tn_p        <- d[[paste0("Tumor_vs_normal_p.adj", suf)]]
    tn_p_num    <- coerce_numeric_safe(tn_p)
    tn_strength <- neg_log10_p(tn_p_num, eps)  # Wilcoxon p, not FDR
    
    # 5–16) Survival blocks for OS, DSS, DFI, PFI
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
      
      # Neutralize strength and χ² whenever direction is 0 (NS)
      strg[dir == 0 | is.na(dir)] <- 0
      chis[dir == 0 | is.na(dir)] <- 0
      
      surv_dir_list[[sname]]      <- dir
      surv_strength_list[[sname]] <- strg
      surv_chisq_list[[sname]]    <- chis
    }
    
    # 17) Microenvironment score
    tme_score <- coerce_numeric_safe(d[[paste0("Microenvironment_score", suf)]])
    tme_score[is.na(tme_score)] <- 0
    
    # 18) Immune direction
    immune_dir <- map_immune_direction(d[[paste0("Immune_classification", suf)]])
    immune_dir[is.na(immune_dir)] <- 0
    
    # Stack in fixed 18D order
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
  
  # Metadata (labeling only)
  meta_cols <- intersect(
    c("Circuitries_id", "CTAB",
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
      "Category", "Signature_count", "Interaction_count"),
    colnames(circuitries)
  )
  
  meta <- as_tibble(circuitries[, meta_cols, drop = FALSE])
  
  if ("Circuitries_id" %in% colnames(circuitries)) {
    rownames(features_sig) <- circuitries$Circuitries_id
    rownames(features_int) <- circuitries$Circuitries_id
    rownames(meta)         <- circuitries$Circuitries_id
  }
  
  if (verbose) {
    message("Finished building latent feature matrices: ",
            nrow(features_sig), " circuitries × ", ncol(features_sig),
            " dimensions per side (18D).")
  }
  
  list(
    features_sig = features_sig,
    features_int = features_int,
    meta         = meta
  )
}

###############################################################################
# 3. PCA embedding from tensor (numerically robust, 18D → 3D)
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
  message("PCA completed. Variance explained by PC1–PC", n_components,
          " = ", round(sum(var_exp) * 100, 2), "%.")
  
  list(
    pca        = pca_obj,
    coords_sig = coords_sig,
    coords_int = coords_int
  )
}

###############################################################################
# 4. Project arbitrary latent vectors into PCA space
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
  
  X_scaled <- scale(X,
                    center = pca_obj$center,
                    scale  = pca_obj$scale)
  scores <- X_scaled %*% pca_obj$rotation
  scores[, seq_len(min(n_components, ncol(scores))), drop = FALSE]
}

###############################################################################
# 5. Build barycenter-centred ± vertices for one latent vector (v-strategy)
###############################################################################
# For v ∈ R^P (P = 18), with coordinates v_j, j = 1..P:
#   step_j = max(|v_j|, ε)   (ε small, e.g. 1e−6)
#   p_{j,+} = v + step_j e_j
#   p_{j,−} = v − step_j e_j
# This yields 2P vertices that define a local axis-aligned shell around v.
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
# 6. Build ONE circuitry polytope: TWO hulls (sig + int) in same space
###############################################################################

build_circuitry_polytope <- function(tensor,
                                     embedding,
                                     index = NULL,
                                     circuitry_id = NULL) {
  
  if (is.null(index) && is.null(circuitry_id)) {
    stop("Provide either 'index' or 'circuitry_id'.")
  }
  
  meta <- tensor$meta
  
  # Resolve by Circuitries_id if provided
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
  
  pca_obj <- embedding$pca
  if (is.null(pca_obj)) {
    stop("embedding$pca is NULL.")
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
  
  # Annotate dimension names (18D)
  verts_sig$dim_name <- latent_dim_names[verts_sig$dim_index]
  verts_int$dim_name <- latent_dim_names[verts_int$dim_index]
  
  vertices_all <- bind_rows(verts_sig, verts_int)
  
  # Barycenters (projected 18D vectors)
  bary_sig_3d <- project_to_pca(pca_obj, matrix(v_sig, nrow = 1), n_components = 3)[1, ]
  bary_int_3d <- project_to_pca(pca_obj, matrix(v_int, nrow = 1), n_components = 3)[1, ]
  names(bary_sig_3d) <- c("PC1", "PC2", "PC3")
  names(bary_int_3d) <- c("PC1", "PC2", "PC3")
  
  # Convex hull for signature side only
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
  
  # Convex hull for interaction side only
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
# 7. Plot circuitry polytope: TWO hulls (sig + int) + 18D legend
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
  
  p <- plot_ly()
  
  # 1) Vertices
  if (show_vertices && nrow(verts_all) > 0L) {
    p <- p %>%
      add_markers(
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
  
  # 2) Hull for sig
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
        add_trace(
          type   = "mesh3d",
          mode = "markers",
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
  
  # 3) Hull for int
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
        add_trace(
          type   = "mesh3d",
          mode = "markers",
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
      add_markers(
        x = bary_sig["PC1"],
        y = bary_sig["PC2"],
        z = bary_sig["PC3"],
        marker   = list(size = 8, symbol = "diamond"),
        name     = "Barycenter_sig",
        text     = "Signature barycenter (18D → 3D)",
        hoverinfo = "text"
      ) %>%
      add_markers(
        x = bary_int["PC1"],
        y = bary_int["PC2"],
        z = bary_int["PC3"],
        marker   = list(size = 8, symbol = "diamond-open"),
        name     = "Barycenter_int",
        text     = "Interaction barycenter (18D → 3D)",
        hoverinfo = "text"
      )
  }
  
  # 5) Text block with the 18 latent dimensions
  dim_legend_text <- paste0(
    "<b>Latent dimensions (18D):</b><br>",
    paste0(seq_along(latent_dim_names), " = ", latent_dim_names, collapse = "<br>")
  )
  
  p <- p %>%
    layout(
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
    saveWidget(as_widget(p), file = safe_file, selfcontained = TRUE)
    message("Polytope HTML saved to: ", normalizePath(safe_file))
  }
  
  p
}

###############################################################################
# 8. Barycenter distance = multidimensional phenotypic discordance
###############################################################################

barycenter_distance <- function(i, tensor, embedding) {
  poly_i <- build_circuitry_polytope(tensor, embedding, index = i)
  sqrt(sum((poly_i$bary_sig_3d - poly_i$bary_int_3d)^2))
}

###############################################################################
# 9. Example usage (assuming `circuitries` already in memory)
###############################################################################
tensor    <- build_geometric_tensor_from_circuitries(circuitries)
embedding <- compute_pca_embedding_from_tensor(tensor, n_components = 3)

# # Example: first circuitry by row index
poly1 <- build_circuitry_polytope(tensor, embedding, index = 17927L)
p1 <- plot_circuitry_polytope(
  poly1,
  outfile_html = paste0(
    "polytope_",
    gsub("[^A-Za-z0-9]+", "_", poly1$circuitry_id),
    "_dual_hulls_18D.html"
  )
)
p1

# Barycenter distance (multidimensional discordance)
# Barycenter distance as multidimensional phenotypic discordance
dist1 <- barycenter_distance(17927L, tensor, embedding)
dist1

# # Or by specific Circuitries_id:
poly_LGG <- build_circuitry_polytope(
  tensor,
  embedding,
  circuitry_id = "LGG-6662 / LGG-5995"
)
plot_circuitry_polytope(poly_LGG)
###############################################################################

## Sanity check for 18D x 2 hulls dimensionallity
poly1 <- build_circuitry_polytope(tensor, embedding, index = 17927L)
ncol(tensor$features_sig)       # should be 18
nrow(poly1$vertices_sig)        # should be 36
nrow(poly1$vertices_int)        # should be 36
nrow(poly1$vertices_all)        # should be 72
unique(poly1$vertices_sig$dim_index)  # 1:18
table(poly1$vertices_sig$dim_index)   # each dim index should appear exactly twice (pos/neg)

## 1) Confirm latent dimensionality
ncol(tensor$features_sig)       # should be 18
ncol(tensor$features_int)       # should be 18

## 2) Build one circuitry polytope
poly1 <- build_circuitry_polytope(tensor, embedding, index = 17927L)

## 3) Check vertex counts
nrow(poly1$vertices_sig)        # should be 36 (2 * 18)
nrow(poly1$vertices_int)        # should be 36
nrow(poly1$vertices_all)        # should be 72

## 4) Check that every latent dimension appears twice (pos/neg) per side
sort(unique(poly1$vertices_sig$dim_index))   # should be 1:18
table(poly1$vertices_sig$dim_index)          # each should be 2
table(poly1$vertices_sig$sign)               # "pos" and "neg" 18 each

sort(unique(poly1$vertices_int$dim_index))   # again 1:18
table(poly1$vertices_int$dim_index)
table(poly1$vertices_int$sign)

# hull volume
# Notes: 
# (1) “Flat” phenotypic behavior; limited regulatory potency; non-specific or single-dimensional signal.”
# (2) Large volume: highly informative, multi-dimensional, and strongly phenotype-coupled.
length(poly1$hull_sig$faces)  # > 0 if hull successfully built
poly1$hull_sig$volume         # finite, non-NA if geometry::convhulln succeeded
poly1$hull_int$volume

###############################################################################
# 8. Barycenter distance = multidimensional phenotypic discordance
###############################################################################
# The barycenter distance is computed in the 3D PCA space as:
#
#   d_bary = || b_sig − b_int ||_2
#
# where:
#   b_sig = (PC1_sig, PC2_sig, PC3_sig)  : barycenter of the signature polytope
#   b_int = (PC1_int, PC2_int, PC3_int)  : barycenter of the interaction polytope
#
# In practice, d_bary behaves as a **multidimensional discordance magnitude**
# between the two 18D latent vectors (signature vs interaction), after global
# PCA embedding.
#
# Heuristic interpretation in PCA units:
#   d_bary < 0.5       → high concordance (very similar latent profiles)
#   0.5–1.5            → mild/moderate discordance
#   > 1.5              → strong discordance
#   > 2.5–3            → extreme discordance
#
# Combined with hull volumes, we obtain a 2D interpretation grid:
#
#   X-axis: barycenter_distance  (discordance magnitude)
#   Y-axis: hull volume          (phenotypic/latent "complexity" of each side)
#
# Let:
#   V_sig = hull_sig$volume   (signature hull volume in PCA space)
#   V_int = hull_int$volume   (interaction hull volume)
#   V_mean = (V_sig + V_int) / 2
#
# Then conceptually:
#
# 1) Small d_bary, small V_mean
#    ---------------------------------
#    - Signature and interaction are close in the latent PCA space.
#    - Both polytopes are compact (small volume).
#    - Interpretation:
#        • Simple, low-dimensional circuitry with highly concordant behavior.
#        • Similar tumor–normal, survival, TME, and immune patterns.
#        • Likely “clean” regulatory module with coherent phenotypic readout.
#
# 2) Small d_bary, large V_mean
#    ---------------------------------
#    - Signature and interaction remain close (concordant barycenters),
#      but each polytope spans a large volume.
#    - Interpretation:
#        • Complex, high-dimensional circuitry with many active latent axes.
#        • Despite internal complexity, the two sides remain phenotypically
#          aligned (high concordance).
#        • Suggests robust, multi-mechanistic coupling between signature and
#          interaction (e.g., multiple survival endpoints + immune/TME axes
#          deforming both hulls in a coordinated way).
#
# 3) Large d_bary, small V_mean
#    ---------------------------------
#    - Signature and interaction barycenters are far apart (discordant),
#      but both hulls are compact.
#    - Interpretation:
#        • Two relatively “simple” but mismatched phenotypic programs.
#        • Each side is low-dimensional and well-defined, yet they encode
#          different directions of effect (e.g., opposite survival/immune
#          directions, contrasting tumor–normal behavior).
#        • This is a clean form of discordance: sharp, low-noise mismatch
#          between two simple modules.
#
# 4) Large d_bary, large V_mean
#    ---------------------------------
#    - Signature and interaction barycenters are far apart.
#    - Both polytopes exhibit high volume (many deformed axes).
#    - Interpretation:
#        • Highly complex, strongly discordant circuitry.
#        • Multiple latent dimensions contribute to the disagreement between
#          the two sides (e.g., divergent survival signals across endpoints,
#          conflicting immune/TME patterns, and opposite tumor–normal shifts).
#        • This regime is a candidate for biologically “chaotic” or
#          context-dependent regulation, where signatures and their interacting
#          partners encode rich but conflicting information about the tumor
#          phenotype.
#
# In summary:
#   - d_bary captures **where** the two sides sit relative to each other in the
#     latent space (concordant vs discordant).
#   - Hull volumes capture **how many latent axes** are effectively contributing
#     to each side (simple vs complex regulatory geometry).
#   - The (d_bary, V_sig, V_int) triplet defines a geometric phenotype of the
#     circuitry: concordant/simple, concordant/complex, discordant/simple,
#     or discordant/complex.
###############################################################################

###############################################################################
# Parallel computation of barycenter distance and convex hull volumes
# for all circuitries (N = nrow(circuitries))
###############################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("parallel", quietly = TRUE))   install.packages("parallel")
  if (!requireNamespace("dplyr", quietly = TRUE))      install.packages("dplyr")
  if (!requireNamespace("geometry", quietly = TRUE))   install.packages("geometry")
})

###############################################################################
# 1. Build latent tensor and PCA embedding (18D → 3D)
###############################################################################

# Assumes the following functions are already defined in the session:
# - build_geometric_tensor_from_circuitries()
# - compute_pca_embedding_from_tensor()
# - build_circuitry_polytope()
# - build_vertices_for_signature()
# - project_to_pca()
# - latent_dim_names (vector)

###############################################################################
# 2. Parallel evaluation of all circuitries (barycenter distance + hull volumes)
###############################################################################

n_rows  <- nrow(circuitries)
indices <- seq_len(n_rows)

# Use 10 cores as requested (adjust if machine has fewer cores)
n_cores <- 10L

cl <- parallel::makeCluster(n_cores)

# Load needed packages and objects on workers
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(geometry)
    library(dplyr)
    library(tibble)
  })
})

parallel::clusterExport(
  cl,
  varlist = c(
    "tensor",
    "embedding",
    "build_circuitry_polytope",
    "build_vertices_for_signature",
    "project_to_pca",
    "latent_dim_names"
  ),
  envir = environment()
)

# Worker function: one circuitry → distance + two hull volumes
compute_geometry_for_index <- function(i) {
  poly_i <- build_circuitry_polytope(tensor, embedding, index = i)
  
  # Barycenter distance (Euclidean in 3D PCA space)
  d_bary <- sqrt(sum((poly_i$bary_sig_3d - poly_i$bary_int_3d)^2))
  
  # Hull volumes (may be NA or non-finite if hull fails)
  v_sig <- poly_i$hull_sig$volume
  v_int <- poly_i$hull_int$volume
  
  if (!is.finite(v_sig)) v_sig <- NA_real_
  if (!is.finite(v_int)) v_int <- NA_real_
  
  list(
    barycenter_distance = d_bary,
    sig_hull_vol        = v_sig,
    int_hull_vol        = v_int
  )
}

# Parallel loop over all circuitries
res_list <- parallel::parLapply(cl, indices, compute_geometry_for_index)

parallel::stopCluster(cl)

###############################################################################
# 3. Unpack results into numeric vectors
###############################################################################

barycenter_distance_vec <- vapply(
  res_list, `[[`, numeric(1), "barycenter_distance"
)
sig_hull_vol_vec <- vapply(
  res_list, `[[`, numeric(1), "sig_hull_vol"
)
int_hull_vol_vec <- vapply(
  res_list, `[[`, numeric(1), "int_hull_vol"
)

# Mean volume per circuitry (phenotypic "complexity" scale)
V_mean_vec <- rowMeans(
  cbind(sig_hull_vol_vec, int_hull_vol_vec),
  na.rm = TRUE
)

# Volume asymmetry ratio: max / min (≈1 → symmetric; >>1 → asymmetric)
eps_vol   <- 1e-12
vol_min   <- pmin(sig_hull_vol_vec, int_hull_vol_vec)
vol_max   <- pmax(sig_hull_vol_vec, int_hull_vol_vec)
vol_ratio_vec <- vol_max / pmax(vol_min, eps_vol)
vol_ratio_vec[!is.finite(vol_ratio_vec)] <- NA_real_

###############################################################################
# 4. Derive biologically interpretable categories from geometric metrics
###############################################################################

# Heuristic breakpoints for barycenter distance (discordance magnitude)
#   < 0.5   → high concordance
#   0.5–1.5 → moderate discordance
#   1.5–2.5 → strong discordance
#   > 2.5   → extreme discordance
distance_breaks <- c(0.5, 1.5, 2.5)

# Data-driven breakpoints for mean hull volume:
#   low_dimensional           → lower tertile
#   intermediate_complexity   → middle tertile
#   high_complexity           → upper tertile
finite_vm <- V_mean_vec[is.finite(V_mean_vec) & V_mean_vec > 0]

if (length(finite_vm) >= 10L) {
  vol_breaks <- as.numeric(stats::quantile(finite_vm, probs = c(1/3, 2/3)))
} else {
  # Fallback if volumes are scarce/degenerate
  vm_med     <- stats::median(finite_vm, na.rm = TRUE)
  vol_breaks <- c(vm_med / 2, vm_med * 2)
}

# Mapping function:
#   Input  : barycenter_distance, vol_ratio, V_mean
#   Output : list(distance_implication, vol_implication)
interpret_geometric_pattern <- function(d_bary, vol_ratio, v_mean) {
  
  # Distance implication: where sig vs int sit in latent space
  distance_implication <- dplyr::case_when(
    !is.finite(d_bary)                  ~ "distance_undefined",
    d_bary < distance_breaks[1]         ~ "high_concordance",
    d_bary < distance_breaks[2]         ~ "moderate_discordance",
    d_bary < distance_breaks[3]         ~ "strong_discordance",
    TRUE                                ~ "extreme_discordance"
  )
  
  # Volume implication: how many latent axes contribute, and how asymmetric
  if (!is.finite(v_mean) || v_mean <= 0) {
    vol_implication <- "flat_or_unresolved_geometry"
  } else {
    # Complexity tier from mean hull volume
    vol_complexity <- dplyr::case_when(
      v_mean <= vol_breaks[1] ~ "low_dimensional_flat",
      v_mean <= vol_breaks[2] ~ "intermediate_complexity",
      TRUE                    ~ "high_complexity_multidimensional"
    )
    
    # Symmetry / asymmetry between sig and int hulls
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

# Vectorised application over all circuitries
imp_list <- Map(
  interpret_geometric_pattern,
  d_bary   = barycenter_distance_vec,
  vol_ratio = vol_ratio_vec,
  v_mean   = V_mean_vec
)

distance_implication_vec <- vapply(
  imp_list, `[[`, character(1), "distance_implication"
)
vol_implication_vec <- vapply(
  imp_list, `[[`, character(1), "vol_implication"
)

###############################################################################
# 5. Attach new variables to `circuitries` in the requested positions
###############################################################################

# Attach all new columns first (temporarily at the end)
circuitries$barycenter_distance <- barycenter_distance_vec
circuitries$distance_implication <- distance_implication_vec

circuitries$sig_hull_vol   <- sig_hull_vol_vec
circuitries$int_hull_vol   <- int_hull_vol_vec
circuitries$vol_ratio      <- vol_ratio_vec
circuitries$vol_implication <- vol_implication_vec

# Column re-ordering:
#   Circuitries_id
#   barycenter_distance
#   distance_implication
#   [all other existing columns, including hull volume metrics]
all_cols <- colnames(circuitries)

# Remove the two distance-related columns temporarily
cols_without_dist <- all_cols[!all_cols %in% c("barycenter_distance",
                                               "distance_implication")]

idx_id <- match("Circuitries_id", cols_without_dist)

new_order <- c(
  cols_without_dist[1:idx_id],
  "barycenter_distance",
  "distance_implication",
  cols_without_dist[(idx_id + 1L):length(cols_without_dist)]
)

circuitries <- circuitries[, new_order]

# Ensure the variables exist before reordering
required_vars <- c("distance_implication",
                   "sig_hull_vol", "int_hull_vol",
                   "vol_ratio", "vol_implication")

missing <- setdiff(required_vars, names(circuitries))
if (length(missing) > 0) {
  stop("Missing variables: ", paste(missing, collapse = ", "))
}

# Current column names
cols <- names(circuitries)

# Position of "distance_implication"
pos_dist_impl <- match("distance_implication", cols)

# New block of columns to insert
vol_block <- c("sig_hull_vol", "int_hull_vol",
               "vol_ratio", "vol_implication")

# Remove these columns from their current positions
cols_remaining <- cols[!cols %in% vol_block]

# Insert the block immediately after "distance_implication"
new_order <- append(cols_remaining, vol_block, after = pos_dist_impl)

# Reorder the data frame
circuitries <- circuitries[, new_order]

###  Sanity check
###  
dist1 <- barycenter_distance(17947L, tensor, embedding)
dist1

# If `circuitries` is already loaded, skip this.
if (!requireNamespace("rio", quietly = TRUE)) install.packages("rio")
library(rio)
circuitries_original <- import("Regulatory_circuitries_02.tsv")

## 1) Identify the circuitry in the CURRENT sorted data frame
row_idx <- 17947L
cid <- circuitries$Circuitries_id[row_idx]

cid
# e.g. "LGG-6708 / LGG-5904" (just an example)

## 2) Extract the distance that is already stored in the data frame
d_df <- circuitries$barycenter_distance[row_idx]

## 3) Recompute the distance by Circuitries_id using the same tensor+embedding
barycenter_distance_by_id <- function(circuitry_id, tensor, embedding) {
  poly_i <- build_circuitry_polytope(
    tensor,
    embedding,
    circuitry_id = circuitry_id
  )
  sqrt(sum((poly_i$bary_sig_3d - poly_i$bary_int_3d)^2))
}

d_recomputed <- barycenter_distance_by_id(cid, tensor, embedding)

d_df
d_recomputed
###############################################################################
# Updated `circuitries` now contains:
##############################################################################
#   - barycenter_distance          (after Circuitries_id)
#   - distance_implication         (after barycenter_distance)
#   - sig_hull_vol, int_hull_vol   (hull volumes for each side)
#   - vol_ratio                    (asymmetry metric)
#   - vol_implication              (biological interpretation of volume pattern)
###############################################################################

# # Example: first circuitry by row index
poly1 <- build_circuitry_polytope(tensor, embedding, index = 17947L)
p1 <- plot_circuitry_polytope(
  poly1,
  outfile_html = paste0(
    "polytope_",
    gsub("[^A-Za-z0-9]+", "_", poly1$circuitry_id),
    "_dual_hulls_18D.html"
  )
)
p1

# Barycenter distance (multidimensional discordance)
# Barycenter distance as multidimensional phenotypic discordance
dist1 <- barycenter_distance(17947, tensor, embedding)
dist1

## 1) Build the polytope object for the selected circuitry
poly_PAAD <- build_circuitry_polytope(
  tensor,
  embedding,
  circuitry_id = "PAAD-3581 / PAAD-7300"
)

plot_circuitry_polytope(poly_PAAD)

poly_LGG <- build_circuitry_polytope(
  tensor,
  embedding,
  circuitry_id = "LGG-6708 / LGG-5904"
)

plot_circuitry_polytope(poly_LGG)

poly_LGG <- build_circuitry_polytope(
  tensor,
  embedding,
  circuitry_id = "LGG-5180 / LGG-2756"
)

plot_circuitry_polytope(poly_LGG)

## 2) Define output filename in the working directory
outfile <- paste0(
  "polytope_",
  gsub("[^A-Za-z0-9_]+", "_", poly_LGG$circuitry_id),
  "_dual_hulls_18D.html"
)

## 3) Generate and save the 3D polytope
p_LGG <- plot_circuitry_polytope(
  poly_LGG,
  outfile_html = outfile     # <-- this saves directly in getwd()
)

## 4) Print confirmation
cat("Polytope saved to:\n", normalizePath(outfile), "\n")

saveRDS(circuitries, file = "Regulatory_circuitries_02_geometric.rds")

rio::export(circuitries, "Regulatory_circuitries_02_geometric.tsv")

##### 
##### 
##### Metabolism and Pathways Overlay Analysis on Geometric Circuitry Classes
##### 
##### 
###############################################################################
# SECTION 1 — Load required packages
###############################################################################
# dplyr   → data manipulation
# tidyr   → reshaping to wide/long formats
# readr   → efficient TSV reading/writing
###############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

###############################################################################
# SECTION 2 — Import precomputed geometric circuitry dataset
###############################################################################
# The file "Regulatory_circuitries_02_geometric.tsv" contains:
#   • geometric variables (barycenter distance, convex-hull implications)
#   • curated metabolic annotations (Pathways, Metabolism)
#   • >100 phenotypic variables (ignored in this analysis)
#
# Only the geometric classification variables and metabolic labels are used.
###############################################################################
circuitries <- read_tsv("Regulatory_circuitries_02_geometric.tsv",
                        show_col_types = FALSE)

###############################################################################
# SECTION 3 — Validate presence of the required analytical variables
###############################################################################
required_vars <- c("Pathways", "Metabolism",
                   "distance_implication", "vol_implication")

missing_vars <- setdiff(required_vars, colnames(circuitries))
if (length(missing_vars) > 0) {
  stop("Missing required variables: ",
       paste(missing_vars, collapse = ", "))
}

###############################################################################
# SECTION 4 — Prepare factor levels with consistent ordering
###############################################################################
# 'distance_implication' defines 4 ordered discordance regimes:
#   1. high_concordance
#   2. moderate_discordance
#   3. strong_discordance
#   4. extreme_discordance
#
# 'vol_implication' has 9 classes describing:
#   • latent complexity (low / intermediate / high)
#   • geometric asymmetry (symmetric / moderately asymmetric / strongly asymmetric)
#
# Ordering ensures reproducibility and deterministic table output.
###############################################################################

distance_levels <- c("high_concordance",
                     "moderate_discordance",
                     "strong_discordance",
                     "extreme_discordance")

vol_levels <- unique(circuitries$vol_implication)

circuitries <- circuitries %>%
  mutate(
    Pathways    = factor(Pathways),
    Metabolism  = factor(Metabolism),
    distance_implication = factor(distance_implication,
                                  levels = distance_levels),
    vol_implication      = factor(vol_implication,
                                  levels = vol_levels)
  )

###############################################################################
# SECTION 5 — Pathways × Distance_implication
###############################################################################
# Computes for each metabolic pathway:
#   • raw counts of circuitries in each discordance regime
#   • within-pathway proportions (n / total within pathway)
#
# This quantifies how each KEGG pathway distributes across the
# geometric discordance classes derived from the latent tensor.
###############################################################################

path_dist_counts <- circuitries %>%
  count(Pathways, distance_implication, name = "n") %>%
  group_by(Pathways) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Wide-format table suitable for manuscript presentation
path_dist_wide <- path_dist_counts %>%
  select(Pathways, distance_implication, prop) %>%
  pivot_wider(
    names_from = distance_implication,
    values_from = prop
  )

###############################################################################
# SECTION 6 — Pathways × Vol_implication
###############################################################################
# Performs the same analysis but across the 9 convex-hull classes,
# characterizing how pathway-level metabolic programs distribute across:
#   • low-dimensional
#   • intermediate-complexity
#   • high-complexity latent geometries
###############################################################################

path_vol_counts <- circuitries %>%
  count(Pathways, vol_implication, name = "n") %>%
  group_by(Pathways) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

path_vol_wide <- path_vol_counts %>%
  select(Pathways, vol_implication, prop) %>%
  pivot_wider(
    names_from = vol_implication,
    values_from = prop
  )

###############################################################################
# SECTION 7 — Metabolism × Distance_implication
###############################################################################
# Aggregates the analysis at the KEGG superfamily level:
#   • Carbohydrate metabolism
#   • Lipid metabolism
#   • Amino acid metabolism
#   • Nucleotide metabolism
#   • Energy metabolism
#   • Cofactors/vitamins metabolism
#   • "Metabolism of other amino acids"
#
# This provides high-level biological structuring across discordance regimes.
###############################################################################

metab_dist_counts <- circuitries %>%
  count(Metabolism, distance_implication, name = "n") %>%
  group_by(Metabolism) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

metab_dist_wide <- metab_dist_counts %>%
  select(Metabolism, distance_implication, prop) %>%
  pivot_wider(
    names_from = distance_implication,
    values_from = prop
  )

###############################################################################
# SECTION 8 — Metabolism × Vol_implication
###############################################################################
# Similar aggregation across convex-hull classes.
# This reveals whether specific metabolic superfamilies preferentially occupy:
#   • low / intermediate / high latent complexity
#   • symmetric vs. asymmetric geometric regimes
###############################################################################

metab_vol_counts <- circuitries %>%
  count(Metabolism, vol_implication, name = "n") %>%
  group_by(Metabolism) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

metab_vol_wide <- metab_vol_counts %>%
  select(Metabolism, vol_implication, prop) %>%
  pivot_wider(
    names_from = vol_implication,
    values_from = prop
  )

###############################################################################
# SECTION 9 — Export analysis tables for manuscript / supplement
###############################################################################
write_tsv(path_dist_counts, "overlay_Pathways_x_distance_counts.tsv")
write_tsv(path_dist_wide,   "overlay_Pathways_x_distance_prop.tsv")

write_tsv(path_vol_counts,  "overlay_Pathways_x_volume_counts.tsv")
write_tsv(path_vol_wide,    "overlay_Pathways_x_volume_prop.tsv")

write_tsv(metab_dist_counts, "overlay_Metabolism_x_distance_counts.tsv")
write_tsv(metab_dist_wide,   "overlay_Metabolism_x_distance_prop.tsv")

write_tsv(metab_vol_counts,  "overlay_Metabolism_x_volume_counts.tsv")
write_tsv(metab_vol_wide,    "overlay_Metabolism_x_volume_prop.tsv")

###############################################################################
# SECTION 10 — Print summaries to console for immediate inspection
###############################################################################
cat("\n===== Pathways × Distance (proportions) =====\n")
print(path_dist_wide)

cat("\n===== Pathways × Volume (proportions) =====\n")
print(path_vol_wide)

cat("\n===== Metabolism × Distance (proportions) =====\n")
print(metab_dist_wide)

cat("\n===== Metabolism × Volume (proportions) =====\n")
print(metab_vol_wide)

####
####
####
####
###############################################################################
# ENTROPY-BASED ORDERING OF METABOLIC PATHWAYS AND SUPERFAMILIES
#
# Rationale:
# ----------
# In the geometric circuitry framework, each metabolic pathway (or metabolic
# superfamily) is represented as a distribution of circuitries across several
# geometric regimes. These regimes include:
#
#   (i)   four distance-based discordance classes, capturing the degree of
#         concordance between signature and interaction barycenters; and
#
#   (ii)  nine volume-based complexity/asymmetry classes, quantifying the
#         latent dimensionality and geometric asymmetry of convex-hull
#         representations in the circuitry space.
#
# Different pathways exhibit different degrees of dispersion across these
# classes. Some pathways are narrowly concentrated in a single geometric regime
# (e.g., always high-concordance; or always high-complexity and symmetric),
# whereas others are broadly distributed across many regimes. This dispersion
# provides a quantitative measure of *geometric pleiotropy*—the extent to which
# a metabolic program participates in multiple architectural phenotypes of the
# circuitry landscape.
#
# To measure this dispersion, we compute Shannon entropy:
#
#       H = − Σ p_i log(p_i)
#
# where p_i is the proportion of circuitries assigned to class i. High entropy
# indicates a broad, multi-regime distribution (high pleiotropy), whereas low
# entropy indicates regime specialization.
#
# Biological Interpretation:
# --------------------------
# Entropy serves as an interpretable statistical proxy for the degree of
# regulatory and phenotypic diversity associated with a metabolic pathway.
#
# • High-entropy pathways tend to participate in heterogeneous geometric
#   configurations. These pathways appear across multiple discordance regimes
#   or complexity classes, suggesting generalized or context-flexible metabolic
#   functions. Such pathways act as *metabolic integrators*, influencing tumor
#   biology in diverse contexts, and frequently overlap with central metabolic
#   axes (e.g., nucleotide metabolism, oxidative phosphorylation).
#
# • Low-entropy pathways exhibit restricted geometric signatures. These are
#   *specialized metabolic programs* that align with a more uniform regulatory
#   behavior. Their circuitries occupy a single geometric regime or a narrow
#   subset, reflecting pathway roles that are condition-specific or tightly
#   constrained by tumor metabolic states.
#
# Why entropy ordering improves visualization:
# --------------------------------------------
# Rather than ordering pathways alphabetically (which has no biological
# meaning), entropy provides a principled, information-theoretic method to
# structure the x-axis. By ordering pathways from highest to lowest entropy:
#
#   • broad, pleiotropic metabolic programs appear on the left;
#   • specialized, narrowly patterned programs appear on the right.
#
# This ordering reveals coherent biological gradients within the stacked
# proportional histograms, transforming each figure into an interpretable map
# of metabolic heterogeneity in geometric circuitry space.
#
# Importantly, entropy is computed separately for:
#
#   • distance classes vs. volume classes
#   • pathways vs. metabolic superfamilies
#
# because dispersion across discordance states is not equivalent to dispersion
# across geometric complexity states. Accordingly, four entropy orderings are
# maintained independently:
#
#   Pathways × Distance
#   Pathways × Volume
#   Metabolism × Distance
#   Metabolism × Volume
#
# preserving the biological specificity of each representation.
#
###############################################################################

###############################################################################
# LOAD REQUIRED PACKAGES
###############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

###############################################################################
# IMPORT DATA
###############################################################################
circuitries <- read_tsv("Regulatory_circuitries_02_geometric.tsv", show_col_types = FALSE)

###############################################################################
# LABEL MAPPINGS (DISTANCE, METABOLISM, PATHWAYS, VOLUME)
###############################################################################

# Distance implication classes
map_distance <- c(
  high_concordance       = "High concordance",
  moderate_discordance   = "Moderate discordance",
  strong_discordance     = "Strong discordance",
  extreme_discordance    = "Extreme discordance"
)

# Metabolism superfamilies
map_metabolism <- c(
  "Carbohydrate metabolism"              = "Carbohydrate",
  "Lipid metabolism"                     = "Lipid",
  "Metabolism of cofactors and vitamins" = "Cofactors/vitamins",
  "Energy metabolism"                    = "Energy",
  "Nucleotide metabolism"                = "Nucleotide",
  "Amino acid metabolism"                = "Amino",
  "Metabolism of other amino acids"      = "Metabolism"
)

# Full approved pathway dictionary
map_pathways <- c(
  "Galactose metabolism"                                = "Galactose",
  "Sphingolipid metabolism"                             = "Sphingolipid",
  "Glycerolipid metabolism"                             = "Glycerolipid",
  "Retinol metabolism"                                  = "Retinol",
  "Glycolysis / Gluconeogenesis"                        = "Glycolysis/Gluconeog",
  "Pyruvate metabolism"                                 = "Pyruvate",
  "One carbon pool by folate"                           = "1C folate pool",
  "Porphyrin metabolism"                                = "Porphyrin",
  "Pantothenate and CoA biosynthesis"                   = "Pantothenate/CoA",
  "Oxidative phosphorylation"                           = "OxPhos",
  "Folate biosynthesis"                                 = "Folate",
  "Purine metabolism"                                   = "Purine",
  "Cysteine and methionine metabolism"                  = "Cys/Met",
  "Butanoate metabolism"                                = "Butanoate",
  "Valine, leucine and isoleucine degradation"          = "BCAA degradation",
  "Pyrimidine metabolism"                               = "Pyrimidine",
  "Steroid biosynthesis"                                = "Steroid",
  "Fatty acid degradation"                              = "FA degradation",
  "Histidine metabolism"                                = "Histidine",
  "Tryptophan metabolism"                               = "Tryptophan",
  "beta-Alanine metabolism"                             = "Beta-alanine",
  "Fatty acid biosynthesis"                             = "FA biosynthesis",
  "Fatty acid elongation"                               = "FA elongation",
  "Lysine degradation"                                  = "Lysine degradation",
  "Biosynthesis of unsaturated fatty acids"             = "Unsat FA biosynth",
  "Ascorbate and aldarate metabolism"                   = "Ascorbate/Aldarate",
  "Propanoate metabolism"                               = "Propanoate",
  "Inositol phosphate metabolism"                       = "Inositol phosphate",
  "Arginine and proline metabolism"                     = "Arg/Pro",
  "Glutathione metabolism"                              = "Glutathione",
  "Selenocompound metabolism"                           = "Selenocompounds",
  "Nicotinate and nicotinamide metabolism"              = "Nicotinate/Niacin",
  "Steroid hormone biosynthesis"                        = "Steroid hormone",
  "Ether lipid metabolism"                              = "Ether lipid",
  "Linoleic acid metabolism"                            = "Linoleic acid",
  "alpha-Linolenic acid metabolism"                     = "Alpha-linolenic acid",
  "Glycerophospholipid metabolism"                      = "Glycerophospholipid",
  "Arginine biosynthesis"                               = "Arg biosynth",
  "Amino sugar and nucleotide sugar metabolism"         = "Amino/nucleotide sugars",
  "Pentose phosphate pathway"                           = "PPP",
  "Thiamine metabolism"                                 = "Thiamine",
  "Riboflavin metabolism"                               = "Riboflavin",
  "Glycine, serine and threonine metabolism"            = "Gly/Ser/Thr",
  "Fructose and mannose metabolism"                     = "Fructose/Mannose",
  "Vitamin B6 metabolism"                               = "VitB6",
  "Lipoic acid metabolism"                              = "Lipoic acid",
  "Pentose and glucuronate interconversions"            = "Pentose/Glucuronate",
  "Arachidonic acid metabolism"                         = "Arachidonic acid",
  "Glyoxylate and dicarboxylate metabolism"             = "Glyoxylate/Dicarboxylate",
  "Sulfur metabolism"                                   = "Sulfur",
  "Citrate cycle (TCA cycle)"                           = "TCA cycle",
  "Phosphonate and phosphinate metabolism"              = "Phosphonate/Phosphinate",
  "Alanine, aspartate and glutamate metabolism"         = "Ala/Asp/Glu",
  "Starch and sucrose metabolism"                       = "Starch/Sucrose",
  "Ubiquinone and other terpenoid-quinone biosynthesis" = "Ubiquinone/Terpenoid",
  "Valine, leucine and isoleucine biosynthesis"         = "BCAA biosynth",
  "Primary bile acid biosynthesis"                      = "Bile acid",
  "Tyrosine metabolism"                                 = "Tyrosine",
  "Phenylalanine, tyrosine and tryptophan biosynthesis" = "Aromatic AA biosynth",
  "Phenylalanine metabolism"                            = "Phenylalanine",
  "Nitrogen metabolism"                                 = "Nitrogen",
  "Taurine and hypotaurine metabolism"                  = "Taurine/Hypotaurine",
  "Biotin metabolism"                                   = "Biotin",
  "D-Amino acid metabolism"                             = "D-Amino acid"
)

# Volume mapping
map_volume <- c(
  "low_dimensional_flat ; symmetric_sig_int"               = "Low – Symmetric",
  "low_dimensional_flat ; moderately_asymmetric_sig_int"   = "Low – Mod asymmetric",
  "low_dimensional_flat ; strongly_asymmetric_sig_int"     = "Low – Strong asymmetric",
  "intermediate_complexity ; symmetric_sig_int"             = "Intermediate – Symmetric",
  "intermediate_complexity ; moderately_asymmetric_sig_int" = "Intermediate – Mod asymmetric",
  "intermediate_complexity ; strongly_asymmetric_sig_int"   = "Intermediate – Strong asymmetric",
  "high_complexity_multidimensional ; symmetric_sig_int"               = "High – Symmetric",
  "high_complexity_multidimensional ; moderately_asymmetric_sig_int"   = "High – Mod asymmetric",
  "high_complexity_multidimensional ; strongly_asymmetric_sig_int"     = "High – Strong asymmetric"
)

###############################################################################
# APPLY ALL LABEL MAPPINGS
###############################################################################
circuitries <- circuitries %>%
  mutate(
    distance_implication = recode(distance_implication, !!!map_distance),
    Metabolism           = recode(Metabolism,           !!!map_metabolism),
    Pathways             = recode(Pathways,             !!!map_pathways),
    vol_implication      = recode(vol_implication,      !!!map_volume)
  )

###############################################################################
# DEFINE FACTOR LEVELS (distance + volume)
###############################################################################
circuitries$distance_implication <- factor(
  circuitries$distance_implication,
  levels = c("High concordance", "Moderate discordance",
             "Strong discordance", "Extreme discordance")
)

circuitries$vol_implication <- factor(
  circuitries$vol_implication,
  levels = c(
    "Low – Symmetric", "Low – Mod asymmetric", "Low – Strong asymmetric",
    "Intermediate – Symmetric", "Intermediate – Mod asymmetric", "Intermediate – Strong asymmetric",
    "High – Symmetric", "High – Mod asymmetric", "High – Strong asymmetric"
  )
)

###############################################################################
# SHANNON ENTROPY FUNCTION
###############################################################################
shannon_entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

###############################################################################
# ENTROPY TABLES (Pathways × Distance / Pathways × Volume)
###############################################################################

path_dist_entropy <- circuitries %>%
  count(Pathways, distance_implication) %>%
  group_by(Pathways) %>%
  mutate(prop = n / sum(n)) %>%
  summarise(entropy = shannon_entropy(prop), .groups = "drop") %>%
  arrange(desc(entropy))

path_vol_entropy <- circuitries %>%
  count(Pathways, vol_implication) %>%
  group_by(Pathways) %>%
  mutate(prop = n / sum(n)) %>%
  summarise(entropy = shannon_entropy(prop), .groups = "drop") %>%
  arrange(desc(entropy))

###############################################################################
# ENTROPY TABLES (Metabolism × Distance / Metabolism × Volume)
###############################################################################

metab_dist_entropy <- circuitries %>%
  count(Metabolism, distance_implication) %>%
  group_by(Metabolism) %>%
  mutate(prop = n / sum(n)) %>%
  summarise(entropy = shannon_entropy(prop), .groups = "drop") %>%
  arrange(desc(entropy))

metab_vol_entropy <- circuitries %>%
  count(Metabolism, vol_implication) %>%
  group_by(Metabolism) %>%
  mutate(prop = n / sum(n)) %>%
  summarise(entropy = shannon_entropy(prop), .groups = "drop") %>%
  arrange(desc(entropy))

###############################################################################
# FOUR ENTROPY-BASED X-AXIS FACTORS
###############################################################################

circuitries$Pathways_dist_order <- factor(
  circuitries$Pathways, levels = path_dist_entropy$Pathways
)

circuitries$Pathways_vol_order <- factor(
  circuitries$Pathways, levels = path_vol_entropy$Pathways
)

circuitries$Metabolism_dist_order <- factor(
  circuitries$Metabolism, levels = metab_dist_entropy$Metabolism
)

circuitries$Metabolism_vol_order <- factor(
  circuitries$Metabolism, levels = metab_vol_entropy$Metabolism
)

###############################################################################
# COLOR PALETTES
###############################################################################
distance_palette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00")

volume_palette <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#D55E00", "#CC79A7", "#F0E442",
  "#0072B2", "#A6761D", "#6A51A3"
)

###############################################################################
# PLOT 1 — Pathways × Distance (entropy-ordered)
###############################################################################
p_path_dist_prop <- ggplot(circuitries,
                           aes(x = Pathways_dist_order,
                               fill = distance_implication)) +
  geom_bar(position = "fill", color = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = distance_palette) +
  labs(x = "Metabolic Pathway", y = "Circuitries (proportion)",
       fill = "Distance Class") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        panel.grid.major.x = element_blank())

tiff("Proportional_Histogram_Pathways_Distance_A4.tiff",
     width = 7016, height = 3508, res = 600, compression = "lzw")
print(p_path_dist_prop)
dev.off()

###############################################################################
# PLOT 2 — Metabolism × Distance (entropy-ordered)
###############################################################################
p_metab_dist_prop <- ggplot(circuitries,
                            aes(x = Metabolism_dist_order,
                                fill = distance_implication)) +
  geom_bar(position = "fill", color = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = distance_palette) +
  labs(x = "Metabolic Superfamily", y = "Circuitries (proportion)",
       fill = "Distance Class") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        panel.grid.major.x = element_blank())

tiff("Proportional_Histogram_Metabolism_Distance_A4.tiff",
     width = 7016, height = 3508, res = 600, compression = "lzw")
print(p_metab_dist_prop)
dev.off()

###############################################################################
# PLOT 3 — Metabolism × Volume (entropy-ordered)
###############################################################################
p_metab_vol_prop <- ggplot(circuitries,
                           aes(x = Metabolism_vol_order,
                               fill = vol_implication)) +
  geom_bar(position = "fill", color = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = volume_palette) +
  labs(x = "Metabolic Superfamily", y = "Circuitries (proportion)",
       fill = "Volume Class") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        panel.grid.major.x = element_blank())

tiff("Proportional_Histogram_Metabolism_Volume_A4.tiff",
     width = 7016, height = 3508, res = 600, compression = "lzw")
print(p_metab_vol_prop)
dev.off()

###############################################################################
# PLOT 4 — Pathways × Volume (entropy-ordered)
###############################################################################
p_path_vol_prop <- ggplot(circuitries,
                          aes(x = Pathways_vol_order,
                              fill = vol_implication)) +
  geom_bar(position = "fill", color = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = volume_palette) +
  labs(x = "Metabolic Pathway", y = "Circuitries (proportion)",
       fill = "Volume Class") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        panel.grid.major.x = element_blank())

tiff("Proportional_Histogram_Pathways_Volume_A4.tiff",
     width = 7016, height = 3508, res = 600, compression = "lzw")
print(p_path_vol_prop)
dev.off()
####
####
####
####

###############################################################################
# ENTROPY ORDERING FOR FINAL_CONCORDANCE_SUMMARY
###############################################################################

# Shannon entropy helper
shannon_entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

# Compute entropy of distance classes within each FCS category
fcs_entropy <- circuitries %>%
  dplyr::count(Final_concordance_summary, distance_implication, name = "n") %>%
  dplyr::group_by(Final_concordance_summary) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::summarise(entropy = shannon_entropy(prop), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(entropy))

print(fcs_entropy)

# Apply entropy-based ordering
circuitries$FCS_entropy_order <- factor(
  circuitries$Final_concordance_summary,
  levels = fcs_entropy$Final_concordance_summary
)

###############################################################################
# COLOR PALETTE (automatically sized)
###############################################################################

fcs_levels <- levels(circuitries$FCS_entropy_order)
k <- length(fcs_levels)

cb_seed <- c(
  "#E69F00", "#56B4E9", "#009E73", "#D55E00",
  "#CC79A7", "#C49A00", "#0072B2", "#A6761D",
  "#6A51A3", "#8C510A", "#4D9221", "#C44E52"
)

fcs_colors  <- grDevices::colorRampPalette(cb_seed)(k)
fcs_palette <- stats::setNames(fcs_colors, fcs_levels)

###############################################################################
# PREPARE DATA FOR PLOT
###############################################################################

fcs_counts <- circuitries %>%
  dplyr::count(FCS_entropy_order, name = "n") %>%
  dplyr::mutate(prop = n / sum(n))

###############################################################################
# PLOT — ENTROPY-ORDERED DISTRIBUTION OF FCS
###############################################################################

p_fcs_entropy <- ggplot2::ggplot(
  fcs_counts,
  ggplot2::aes(x = FCS_entropy_order, y = prop, fill = FCS_entropy_order)
) +
  ggplot2::geom_col(color = "grey20", linewidth = 0.5) +
  ggplot2::scale_fill_manual(values = fcs_palette, drop = FALSE) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.04))
  ) +
  ggplot2::scale_x_discrete(
    expand = ggplot2::expansion(mult = c(0.04, 0.01))
  ) +
  ggplot2::labs(
    x = "Final Concordance Summary (Entropy-ordered)",
    y = "Proportion of Circuitries"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 60,        # increased rotation
      hjust = 1,
      vjust = 1,
      size = 10,
      margin = ggplot2::margin(t = 4)
    ),
    axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 8)
    ),
    axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 12)
    ),
    plot.margin = ggplot2::margin(
      t = 16,
      r = 18,
      b = 32,
      l = 65,
      unit = "pt"
    ),
    legend.position = "none"
  )

###############################################################################
# EXPORT — TIFF, LEGAL LANDSCAPE
###############################################################################

grDevices::tiff(
  filename = "Figure_S3_Distribution_FCS_EntropyOrdered_Legal_600dpi.tiff",
  width = 14,
  height = 8.5,
  units = "in",
  res = 600,
  compression = "lzw",
  bg = "white"
)

print(p_fcs_entropy)
grDevices::dev.off()

###############################################################################
# Load data
###############################################################################
library(dplyr)
library(readr)

circuitries <- read_tsv(
  "Regulatory_circuitries_02_geometric.tsv",
  show_col_types = FALSE
)

###############################################################################
# Clean and ensure numeric structure
###############################################################################
circuitries <- circuitries %>%
  mutate(
    barycenter_distance = as.numeric(barycenter_distance),
    sig_hull_vol        = as.numeric(sig_hull_vol),
    int_hull_vol        = as.numeric(int_hull_vol)
  ) %>%
  mutate(
    mean_volume = rowMeans(cbind(sig_hull_vol, int_hull_vol), na.rm = TRUE)
  )

###############################################################################
# Filter EXTREME discordance class only
###############################################################################
extreme_df <- circuitries %>%
  filter(distance_implication == "extreme_discordance") %>%
  filter(is.finite(barycenter_distance)) %>%
  arrange(desc(barycenter_distance))

###############################################################################
# Identify TOP-N highly discordant circuitries
# Criteria:
#   - Largest barycenter distance
#   - Large geometric mean hull volume
#   - Prefer "High – Mod asymmetric" or "High – Strong asymmetric"
###############################################################################

top_discordant <- extreme_df %>%
  mutate(
    volume_rank = dense_rank(desc(mean_volume)),
    distance_rank = dense_rank(desc(barycenter_distance))
  ) %>%
  mutate(
    composite_score = distance_rank + volume_rank
  ) %>%
  arrange(composite_score) %>%
  slice_head(n = 10)               # <-- Select top 10 discordant circuitries

###############################################################################
# Extract key descriptors for text reporting
###############################################################################

discordant_summary <- top_discordant %>%
  transmute(
    Circuitry_ID      = Circuitries_id,
    Signature         = Nomenclature_sig,
    Interaction       = Nomenclature_int,
    BarycenterDist    = round(barycenter_distance, 4),
    SigVol            = round(sig_hull_vol, 4),
    IntVol            = round(int_hull_vol, 4),
    MeanVol           = round(mean_volume, 4),
    VolumePattern     = vol_implication,
    Metabolism        = Metabolism,
    Pathway           = Pathways
  )

print(discordant_summary)

library(dplyr)

df0 <- circuitries %>%
  mutate(
    vol_ratio_num = as.numeric(vol_ratio),
    sig_vol = as.numeric(sig_hull_vol),
    int_vol = as.numeric(int_hull_vol)
  )

nice_discordant <- df0 %>%
  filter(
    distance_implication %in% c("strong_discordance", "extreme_discordance"),
    barycenter_distance > 2,
    is.finite(vol_ratio_num),
    vol_ratio_num >= 1, vol_ratio_num <= 5,
    is.finite(sig_vol), sig_vol > 0,
    is.finite(int_vol), int_vol > 0
  ) %>%
  mutate(
    mean_vol = (sig_vol + int_vol)/2,
    score = barycenter_distance * mean_vol
  ) %>%
  arrange(desc(score)) %>%
  slice_head(n = 10)

nice_discordant

###############################################################################
# Batch generation of polytope plots for the 10 selected discordant circuitries
###############################################################################

# Ensure nice_discordant has a column named Circuitries_id
stopifnot("Circuitries_id" %in% colnames(nice_discordant))

# Extract Circuitries_id vector
selected_ids <- nice_discordant$Circuitries_id

# Directory to save HTML files (optional)
outdir <- "Discordant_Polytopes"
if (!dir.exists(outdir)) dir.create(outdir)

# Loop over each selected circuitry
for (cid in selected_ids) {
  
  message("Building polytope for circuitry: ", cid)
  
  # Build polytope object
  poly_obj <- build_circuitry_polytope(
    tensor,
    embedding,
    circuitry_id = cid
  )
  
  # Safe filename mapping
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", cid)
  outfile <- file.path(outdir, paste0("polytope_", safe_name, "_dual_hulls_18D.html"))
  
  # Plot and save
  p <- plot_circuitry_polytope(
    poly_obj,
    outfile_html = outfile
  )
  
  # Optionally: also print to viewer
  print(p)
}

message("Completed polytope generation for all selected discordant circuitries.")

###############################################################################
# Select aesthetically ideal discordant circuitries:
#   - strong OR extreme discordance
#   - both hull volumes >= 5   (adjust as desired)
#   - volume ratio between 1 and 2  (balanced geometry)
###############################################################################

ideal_discordant <- circuitries %>%
  filter(
    distance_implication %in% c("strong_discordance", "extreme_discordance"),
    is.finite(sig_hull_vol),
    is.finite(int_hull_vol),
    sig_hull_vol >= 5,
    int_hull_vol >= 5,
    vol_ratio >= 1,
    vol_ratio <= 2
  ) %>%
  arrange(desc(barycenter_distance))   # show strongest discordance first

# Inspect the selected circuitries
ideal_discordant %>%
  select(Circuitries_id,
         barycenter_distance,
         sig_hull_vol,
         int_hull_vol,
         vol_ratio,
         distance_implication) %>%
  head(20)

selected_ids <- ideal_discordant$Circuitries_id

outdir <- "Ideal_Discordant_Polytopes"
if (!dir.exists(outdir)) dir.create(outdir)

for (cid in selected_ids) {
  
  message("Plotting ideal discordant circuitry: ", cid)
  
  poly_obj <- build_circuitry_polytope(tensor, embedding, circuitry_id = cid)
  
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", cid)
  outfile <- file.path(outdir, paste0("polytope_", safe_name, "_dual_hulls_18D.html"))
  
  p <- plot_circuitry_polytope(poly_obj, outfile_html = outfile)
  print(p)
}

library (rio)
rio::export(ideal_discordant, "ideal discordant circuitries.tsv")

df0[df0$Circuitries_id == "LGG-3207 / LGG-2427", ]

df0_01 <- df0 %>% 
  filter(Circuitries_id == "LGG-3207 / LGG-2427")

rio::export(df0_01, "df0_01.tsv")

# ------------------------------------------------------------------------------
# Extract all variables required to write the figure legend
# ------------------------------------------------------------------------------

# Choose circuitry ID
target_id <- "LGG-3207 / LGG-2427"

# Subset the row
cir_3207_2427 <- df0_01[df0_01$Circuitries_id == target_id, ]

# Verify we have exactly one row
if (nrow(cir_3207_2427) != 1) {
  stop("Circuitry ID not found or duplicated.")
}

library(dplyr)
library(tidyr)

## transpose data
cir_tidy_3207_2427 <- cir_3207_2427 %>%
  mutate(across(everything(), ~ as.character(.))) %>%   # force all columns to character
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Value"
  )

cir_tidy_3207_2427

rio::export(cir_tidy_3207_2427, "cir_tidy_3207_2427.tsv")

###############################################################################
# EXTRACT RELEVANT VARIABLES FOR THE CIRCUITRY: LGG-3207 / LGG-2427
###############################################################################

# 1. Target circuitry identifier
target_id <- "LGG-3207 / LGG-2427"

# 2. Variables to extract (same structure as cir_relevant_variables.tsv)
vars_to_extract <- c(
  "Circuitries_id",
  "Nomenclature_sig",
  "Signatures",
  "Nomenclature_int",
  "Interaction",
  "barycenter_distance",
  "sig_hull_vol",
  "int_hull_vol",
  "vol_ratio",
  "distance_implication",
  "vol_implication",
  "Omic_layer_sig",
  "Phenotypic_layer_sig",
  "Omic_layer_int",
  "Phenotypic_layer_int",
  "Correlation_rho_sig",
  "Correlation_rho_int",
  "Tumor_vs_normal_sig",
  "Tumor_vs_normal_int",
  "Combined_outcome_HRC_sig",
  "Combined_outcome_HRC_int",
  "Microenvironment_classification_sig",
  "Microenvironment_classification_int",
  "Immune_classification_sig",
  "Immune_classification_int",
  "Phenotypic_concordance",
  "Survival_concordance_aggregated",
  "Immune_concordance",
  "Final_concordance_summary",
  "CTAB",
  "Metabolism",
  "Pathways"
)
# 3. Extract, convert all values to character (required for tidy pivot), pivot to two-column format
cir_3207_2427_relevant_variables <- circuitries %>%
  dplyr::filter(Circuitries_id == target_id) %>%
  dplyr::select(dplyr::all_of(vars_to_extract)) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%  # ensure homogeneous type
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# 4. Print result
cir_3207_2427_relevant_variables

# 5. (Optional) Save to TSV
readr::write_tsv(cir_3207_2427_relevant_variables, "LGG_3207_LGG_2427_cir_relevant_variables.tsv")
###############################################################################

###############################################################################
# EXTRACT RELEVANT VARIABLES FOR THE CIRCUITRY: LGG-6708 / LGG-5904
###############################################################################

# 1. Target circuitry identifier
target_id <- "LGG-6708 / LGG-5904"

# 2. Variables to extract (same structure as cir_relevant_variables.tsv)
vars_to_extract <- c(
  "Circuitries_id",
  "Nomenclature_sig",
  "Signatures",
  "Nomenclature_int",
  "Interaction",
  "barycenter_distance",
  "sig_hull_vol",
  "int_hull_vol",
  "vol_ratio",
  "distance_implication",
  "vol_implication",
  "Omic_layer_sig",
  "Phenotypic_layer_sig",
  "Omic_layer_int",
  "Phenotypic_layer_int",
  "Correlation_rho_sig",
  "Correlation_rho_int",
  "Tumor_vs_normal_sig",
  "Tumor_vs_normal_int",
  "Combined_outcome_HRC_sig",
  "Combined_outcome_HRC_int",
  "Microenvironment_classification_sig",
  "Microenvironment_classification_int",
  "Immune_classification_sig",
  "Immune_classification_int",
  "Phenotypic_concordance",
  "Survival_concordance_aggregated",
  "Immune_concordance",
  "Final_concordance_summary",
  "CTAB",
  "Metabolism",
  "Pathways"
)
# 3. Extract, convert all values to character (required for tidy pivot), pivot to two-column format
cir_6708_5904_relevant_variables <- circuitries %>%
  dplyr::filter(Circuitries_id == target_id) %>%
  dplyr::select(dplyr::all_of(vars_to_extract)) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%  # ensure homogeneous type
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# 4. Print result
cir_6708_5904_relevant_variables

# 5. (Optional) Save to TSV
readr::write_tsv(cir_6708_5904_relevant_variables, "LGG-6708_LGG-5904_cir_relevant_variables.tsv")
###############################################################################

###############################################################################
# OVERRIDE PLOT AESTHETICS FOR HULL COLORS, BARYCENTER SYMBOLS, AND VERTEX SYMBOLS
# + RESTORED 18-DIMENSION LEGEND ANNOTATION (from canonical function)
###############################################################################

plot_circuitry_polytope_custom <- function(poly_obj,
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
  
  ###########################################################################
  # Correct 18D dimension labels — authoritative vector
  ###########################################################################
  latent_dim_names <- poly_obj$latent_dim_names
  
  if (is.null(latent_dim_names) || length(latent_dim_names) != 18) {
    latent_dim_names <- c(
      "rho","rho_strength",
      "TN_dir","TN_strength",
      "OS_dir","OS_strength","OS_lr_chisq",
      "DSS_dir","DSS_strength","DSS_lr_chisq",
      "DFI_dir","DFI_strength","DFI_lr_chisq",
      "PFI_dir","PFI_strength","PFI_lr_chisq",
      "TME_score","Immune_dir"
    )
  }
  
  if (is.null(title)) {
    title <- paste0("Circuitry polytope (18D → 3D): ", poly_obj$circuitry_id)
  }
  
  p <- plot_ly()
  
  ###############################################################################
  # 1 — Vertices (asterisk = positive; circle-open = negative)
  ###############################################################################
  if (show_vertices && nrow(verts_all) > 0L) {
    
    verts_all <- verts_all %>%
      mutate(symbol = ifelse(sign == "pos", "asterisk-open", "circle-open"))
    
    p <- p %>%
      add_markers(
        data   = verts_all,
        x      = ~PC1, y = ~PC2, z = ~PC3,
        symbol = ~symbol,
        symbols = c("asterisk-open", "circle-open"),
        marker = list(size = 5, line = list(width = 1)),
        color  = ~side,
        colors = c("sig" = "darkgreen", "int" = "darkred"),
        text   = ~paste("Side:", side,
                        "<br>Dimension:", dim_name,
                        "<br>Sign:", sign),
        hoverinfo = "text",
        name      = "Vertices"
      )
  }
  
  ###############################################################################
  # 2 — Signature hull (GREEN)
  ###############################################################################
  if (show_hull && !is.null(hull_sig$faces)) {
    
    simplices <- hull_sig$faces
    coords <- as.matrix(verts_sig[, c("PC1", "PC2", "PC3")])
    
    tri <- as.vector(t(simplices)) - 1L
    n_tri <- nrow(simplices)
    
    p <- p %>% add_trace(
      type   = "mesh3d",
      mode = "markers",
      x      = coords[,1], y = coords[,2], z = coords[,3],
      i      = tri[seq(1, 3*n_tri, 3)],
      j      = tri[seq(2, 3*n_tri, 3)],
      k      = tri[seq(3, 3*n_tri, 3)],
      opacity = 0.35,
      color   = "green",
      name    = "Signature hull"
    )
  }
  
  ###############################################################################
  # 3 — Interaction hull (ORANGE)
  ###############################################################################
  if (show_hull && !is.null(hull_int$faces)) {
    
    simplices <- hull_int$faces
    coords <- as.matrix(verts_int[, c("PC1", "PC2", "PC3")])
    
    tri <- as.vector(t(simplices)) - 1L
    n_tri <- nrow(simplices)
    
    p <- p %>% add_trace(
      type   = "mesh3d",
      mode = "markers",
      x      = coords[,1], y = coords[,2], z = coords[,3],
      i      = tri[seq(1, 3*n_tri, 3)],
      j      = tri[seq(2, 3*n_tri, 3)],
      k      = tri[seq(3, 3*n_tri, 3)],
      opacity = 0.35,
      color   = "orange",
      name    = "Interaction hull"
    )
  }
  
  ###############################################################################
  # 4 — Barycenters (dark green & dark red diamonds)
  ###############################################################################
  if (show_barycenters) {
    p <- p %>%
      add_markers(
        x = bary_sig["PC1"], y = bary_sig["PC2"], z = bary_sig["PC3"],
        marker = list(size = 10, symbol = "diamond", color = "darkgreen"),
        name = "Signature barycenter"
      ) %>%
      add_markers(
        x = bary_int["PC1"], y = bary_int["PC2"], z = bary_int["PC3"],
        marker = list(size = 10, symbol = "diamond", color = "darkred"),
        name = "Interaction barycenter"
      )
  }
  
  ###############################################################################
  # 5 — Base layout
  ###############################################################################
  p <- p %>% layout(
    title = list(text = title),
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )
  
  ###############################################################################
  # 6 — RESTORED 18D LEGEND
  ###############################################################################
  
  dim_legend_text <- paste0(
    "<b>Latent dimensions (18D):</b><br>",
    paste0(seq_along(latent_dim_names), " = ", latent_dim_names, collapse = "<br>")
  )
  
  p <- p %>% layout(
    annotations = list(
      list(
        x = 0.02, y = 0.98,
        xref = "paper", yref = "paper",
        showarrow = FALSE,
        align = "left",
        text = dim_legend_text,
        bordercolor = "black",
        borderwidth = 0.5,
        bgcolor = "rgba(255,255,255,0.75)"
      )
    )
  )
  
  ###############################################################################
  # 7 — Save HTML if requested
  ###############################################################################
  if (!is.null(outfile_html)) {
    htmlwidgets::saveWidget(as_widget(p), outfile_html, selfcontained = TRUE)
  }
  
  p
}

###############################################################################
# USAGE EXAMPLES
###############################################################################

poly_LGG <- build_circuitry_polytope(
  tensor, embedding,
  circuitry_id = "LGG-3207 / LGG-2427"
)

plot_circuitry_polytope_custom(
  poly_LGG,
  outfile_html = "LGG_3207_2427_custom.html"
)

poly_LGG_1 <- build_circuitry_polytope(
  tensor, embedding,
  circuitry_id = "LGG-6708 / LGG-5904"
)

plot_circuitry_polytope_custom(
  poly_LGG_1,
  outfile_html = "LGG_6708_5904_custom.html"
)

#### Aesthetics examples
###############################################################################
# OVERRIDE PLOT AESTHETICS FOR HULL COLORS, BARYCENTER SYMBOLS, AND VERTEX SYMBOLS
# + RESTORED 18-DIMENSION LEGEND ANNOTATION (from canonical function)
###############################################################################

plot_circuitry_polytope_custom <- function(poly_obj,
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
  
  ###########################################################################
  # Correct 18D dimension labels — authoritative vector
  ###########################################################################
  latent_dim_names <- poly_obj$latent_dim_names
  
  if (is.null(latent_dim_names) || length(latent_dim_names) != 18) {
    latent_dim_names <- c(
      "rho","rho_strength",
      "TN_dir","TN_strength",
      "OS_dir","OS_strength","OS_lr_chisq",
      "DSS_dir","DSS_strength","DSS_lr_chisq",
      "DFI_dir","DFI_strength","DFI_lr_chisq",
      "PFI_dir","PFI_strength","PFI_lr_chisq",
      "TME_score","Immune_dir"
    )
  }
  
  if (is.null(title)) {
    title <- paste0("Circuitry polytope (18D → 3D): ", poly_obj$circuitry_id)
  }
  
  p <- plot_ly()
  
  ###############################################################################
  # 1 — Vertices (asterisk = positive; circle-open = negative)
  ###############################################################################
  if (show_vertices && nrow(verts_all) > 0L) {
    
    verts_all <- verts_all %>%
      mutate(symbol = ifelse(sign == "pos", "asterisk-open", "circle-open"))
    
    p <- p %>%
      add_markers(
        data   = verts_all,
        x      = ~PC1, y = ~PC2, z = ~PC3,
        symbol = ~symbol,
        symbols = c("asterisk-open", "circle-open"),
        marker = list(size = 5, line = list(width = 1)),
        color  = ~side,
        colors = c("sig" = "darkblue", "int" = "darkorange"),
        text   = ~paste("Side:", side,
                        "<br>Dimension:", dim_name,
                        "<br>Sign:", sign),
        hoverinfo = "text",
        name      = "Vertices"
      )
  }
  ###############################################################################
  # 2 — Signature hull (LIGHT GREEN)
  ###############################################################################
  if (show_hull && !is.null(hull_sig$faces)) {
    
    simplices <- hull_sig$faces
    coords <- as.matrix(verts_sig[, c("PC1", "PC2", "PC3")])
    
    tri <- as.vector(t(simplices)) - 1L
    n_tri <- nrow(simplices)
    
    # force flat light-green
    vert_col_sig <- rep("#2CA9BC", nrow(coords))   # light green
    
    p <- p %>% add_trace(
      type   = "mesh3d",
      mode = "markers",
      x      = coords[,1], 
      y      = coords[,2], 
      z      = coords[,3],
      i      = tri[seq(1, 3*n_tri, 3)],
      j      = tri[seq(2, 3*n_tri, 3)],
      k      = tri[seq(3, 3*n_tri, 3)],
      opacity     = 0.35,
      vertexcolor = vert_col_sig,
      flatshading = TRUE,
      name        = "Signature hull"
    )
  }
  
  ###############################################################################
  # 3 — Interaction hull (LIGHT ORANGE)
  ###############################################################################
  if (show_hull && !is.null(hull_int$faces)) {
    
    simplices <- hull_int$faces
    coords <- as.matrix(verts_int[, c("PC1", "PC2", "PC3")])
    
    tri <- as.vector(t(simplices)) - 1L
    n_tri <- nrow(simplices)
    
    # force flat light-orange
    vert_col_int <- rep("#F2A65A", nrow(coords))   # light orange
    
    p <- p %>% add_trace(
      type   = "mesh3d",
      mode = "markers",
      x      = coords[,1], 
      y      = coords[,2], 
      z      = coords[,3],
      i      = tri[seq(1, 3*n_tri, 3)],
      j      = tri[seq(2, 3*n_tri, 3)],
      k      = tri[seq(3, 3*n_tri, 3)],
      opacity     = 0.35,
      vertexcolor = vert_col_int,
      flatshading = TRUE,
      name        = "Interaction hull"
    )
  }
  
  ###############################################################################
  # 4 — Barycenters (dark green & dark red diamonds)
  ################################################################################
  if (show_barycenters) {
    p <- p %>%
      add_markers(
        x = bary_sig["PC1"], y = bary_sig["PC2"], z = bary_sig["PC3"],
        marker = list(size = 10, symbol = "diamond", color = "darkblue"),
        name = "Signature barycenter"
      ) %>%
      add_markers(
        x = bary_int["PC1"], y = bary_int["PC2"], z = bary_int["PC3"],
        marker = list(size = 10, symbol = "diamond", color = "darkorange"),
        name = "Interaction barycenter"
      )
  }
  
  ###############################################################################
  # 5 — Base layout
  ###############################################################################
  p <- p %>% layout(
    title = list(text = title),
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )
  
  ###############################################################################
  # 6 — RESTORED 18D LEGEND
  ###############################################################################
  
  dim_legend_text <- paste0(
    "<b>Latent dimensions (18D):</b><br>",
    paste0(seq_along(latent_dim_names), " = ", latent_dim_names, collapse = "<br>")
  )
  
  p <- p %>% layout(
    annotations = list(
      list(
        x = 0.02, y = 0.98,
        xref = "paper", yref = "paper",
        showarrow = FALSE,
        align = "left",
        text = dim_legend_text,
        bordercolor = "black",
        borderwidth = 0.5,
        bgcolor = "rgba(255,255,255,0.75)"
      )
    )
  )
  
  ###############################################################################
  # 7 — Save HTML if requested
  ###############################################################################
  if (!is.null(outfile_html)) {
    htmlwidgets::saveWidget(as_widget(p), outfile_html, selfcontained = TRUE)
  }
  
  p
}

###############################################################################
# USAGE EXAMPLES
###############################################################################

#### Figure 4A. Geometric dual-polytope representation of a high-concordance metabolic regulatory circuitry (LGG-6708 / LGG-5904: SNHG12 → KMT2B). 
plot_circuitry_polytope_custom(
  poly_LGG_1,
  outfile_html = "LGG_6708_5904_custom.html"
)

#### Figure 4A. Geometric dual-polytope representation of a high-concordance metabolic regulatory circuitry (LGG-6708 / LGG-5904: SNHG12 → KMT2B). 
poly_LGG_1 <- build_circuitry_polytope(
  tensor, embedding,
  circuitry_id = "LGG-6708 / LGG-5904"
)

#### Figure 3. Geometric dual-polytope representation of a highly discordant metabolic regulatory circuitry (LGG-3207 / LGG-2427: miR-10a-5p → CYP8B1). 
poly_LGG <- build_circuitry_polytope(
  tensor, embedding,
  circuitry_id = "LGG-3207 / LGG-2427"
)

#### Figure 3. Geometric dual-polytope representation of a highly discordant metabolic regulatory circuitry (LGG-3207 / LGG-2427: miR-10a-5p → CYP8B1). 
plot_circuitry_polytope_custom(
  poly_LGG,
  outfile_html = "LGG_3207_2427_custom.html"
)

###############################################################################
# Figure 2. Conceptual geometric representations of omic signature organization.
# Conceptual geometric figures for Section 1.4 (Figure 2A–F)
# Synthetic data only — NO biological content
# Output: individual plotly figs + one combined multipanel HTML
###############################################################################

set.seed(42)

library(plotly)
library(geometry)

# -----------------------------
# 1) Synthetic point clouds
# -----------------------------
generate_point_cloud <- function(type = c("compact",
                                          "elongated",
                                          "single_axis",
                                          "multi_axis",
                                          "irregular"),
                                 n = 30,
                                 scale = 1) {
  type <- match.arg(type)
  
  if (type == "compact") {
    # isotropic / coherent
    X <- matrix(rnorm(n * 3, sd = 0.35), ncol = 3)
    
  } else if (type == "elongated") {
    # elongated along one axis (anisotropic)
    X <- cbind(
      rnorm(n, sd = 1.20),
      rnorm(n, sd = 0.25),
      rnorm(n, sd = 0.25)
    )
    
  } else if (type == "single_axis") {
    # very strong single-axis dominance
    X <- cbind(
      rnorm(n, sd = 1.50),
      rnorm(n, sd = 0.12),
      rnorm(n, sd = 0.12)
    )
    
  } else if (type == "multi_axis") {
    # two-axis structure (directional multi-axis)
    X <- cbind(
      rnorm(n, sd = 1.00),
      rnorm(n, sd = 0.95),
      rnorm(n, sd = 0.15)
    )
    
  } else if (type == "irregular") {
    # heterogeneous / modular: mixture of clusters
    n1 <- floor(n / 3)
    n2 <- floor(n / 3)
    n3 <- n - n1 - n2
    
    X1 <- matrix(rnorm(n1 * 3, mean = -1.2, sd = 0.25), ncol = 3)
    X2 <- matrix(rnorm(n2 * 3, mean =  0.0, sd = 0.30), ncol = 3)
    X3 <- matrix(rnorm(n3 * 3, mean =  1.2, sd = 0.25), ncol = 3)
    X  <- rbind(X1, X2, X3)
  }
  
  X * scale
}

# -----------------------------
# 2) Convex hull computation
# -----------------------------
compute_hull <- function(X) {
  # options="FA" returns a list with $hull (facets), $area, $vol, etc.
  ch <- convhulln(X, options = "FA")
  
  list(
    points = X,
    facets = ch$hull,   # matrix with 3 columns (triangular facets)
    area   = ch$area,
    vol    = ch$vol
  )
}

# -----------------------------
# 3) Plotly 3D scatter + hull mesh
# -----------------------------
plot_hull_3d <- function(hull_obj,
                         title = NULL,
                         point_color = "gray40",
                         hull_color  = "rgba(100,100,200,0.25)") {
  
  X <- hull_obj$points
  F <- hull_obj$facets
  
  stopifnot(is.matrix(X), ncol(X) == 3)
  stopifnot(is.matrix(F), ncol(F) == 3)
  
  plot_ly() %>%
    add_trace(
      x = X[,1], y = X[,2], z = X[,3],
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 3, color = point_color),
      name = "Components",
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "mesh3d",
      mode = "markers",
      x = X[,1], y = X[,2], z = X[,3],
      i = F[,1] - 1, j = F[,2] - 1, k = F[,3] - 1,  # plotly uses 0-based indexing
      opacity = 0.25,
      color = hull_color,
      name = "Convex hull",
      showlegend = FALSE
    ) %>%
    layout(
      title = list(text = title, x = 0.02, xanchor = "left"),
      margin = list(l = 0, r = 0, b = 0, t = 40),
      scene = list(
        xaxis = list(title = "Latent axis 1"),
        yaxis = list(title = "Latent axis 2"),
        zaxis = list(title = "Latent axis 3")
      )
    )
}

# -----------------------------
# 4) Build all conceptual panels (A–F)
# -----------------------------
# Figure X panels for Section 1.4 (conceptual; simulated only)
X_A <- generate_point_cloud("compact",     n = 30)
X_B <- generate_point_cloud("compact",     n = 30)
X_C <- generate_point_cloud("elongated",   n = 30)
X_D <- generate_point_cloud("single_axis", n = 30)
X_E <- generate_point_cloud("multi_axis",  n = 30)
X_F <- generate_point_cloud("irregular",   n = 30)

h_A <- compute_hull(X_A)
h_B <- compute_hull(X_B)
h_C <- compute_hull(X_C)
h_D <- compute_hull(X_D)
h_E <- compute_hull(X_E)
h_F <- compute_hull(X_F)

fig_A <- plot_hull_3d(h_A, title = "(A) Conceptual projection of a multidimensional omic signature")
fig_B <- plot_hull_3d(h_B, title = "(B) Compact and isotropic geometric organization")
fig_C <- plot_hull_3d(h_C, title = "(C) Elongated and anisotropic geometric organization")
fig_D <- plot_hull_3d(h_D, title = "(D) Coherent single-axis regulatory organization")
fig_E <- plot_hull_3d(h_E, title = "(E) Directional multi-axis regulatory organization")
fig_F <- plot_hull_3d(h_F, title = "(F) Heterogeneous and modular organizational structure")

# Print any single panel interactively, e.g.:
fig_A
fig_B
fig_C
fig_D
fig_E
fig_F
# 
# # -----------------------------
# # 5) Combined multipanel layout (2×3) as one HTML
# # -----------------------------
library(plotly)
library(htmlwidgets)

make_multiscene_3d <- function(fig_list, nrows = 2, ncols = 3,
                               margin_x = 0.03, margin_y = 0.06) {
  
  stopifnot(length(fig_list) == nrows * ncols)
  
  x_breaks <- seq(0, 1, length.out = ncols + 1)
  y_breaks <- seq(0, 1, length.out = nrows + 1)
  
  x_domains <- lapply(seq_len(ncols), function(j) {
    c(x_breaks[j] + margin_x/2, x_breaks[j+1] - margin_x/2)
  })
  
  y_domains <- lapply(seq_len(nrows), function(i) {
    top_row_index <- nrows - i + 1
    c(y_breaks[top_row_index] + margin_y/2, y_breaks[top_row_index+1] - margin_y/2)
  })
  
  out <- plot_ly()
  out$x$data <- list()
  
  for (k in seq_along(fig_list)) {
    
    pb <- plotly_build(fig_list[[k]])
    scene_id <- if (k == 1) "scene" else paste0("scene", k)
    
    for (tr in pb$x$data) {
      tr$scene <- scene_id
      out$x$data <- c(out$x$data, list(tr))
    }
  }
  
  lay <- list(
    margin = list(l = 10, r = 10, b = 10, t = 40)
  )
  
  for (k in seq_along(fig_list)) {
    
    scene_id <- if (k == 1) "scene" else paste0("scene", k)
    
    row_i <- ceiling(k / ncols)
    col_j <- k - (row_i - 1) * ncols
    
    pb <- plotly_build(fig_list[[k]])
    this_title <- ""
    if (!is.null(pb$x$layout$title) && !is.null(pb$x$layout$title$text)) {
      this_title <- pb$x$layout$title$text
    }
    
    lay[[scene_id]] <- list(
      domain = list(x = x_domains[[col_j]], y = y_domains[[row_i]]),
      xaxis  = list(title = "Latent axis 1"),
      yaxis  = list(title = "Latent axis 2"),
      zaxis  = list(title = "Latent axis 3")
    )
    
    if (nzchar(this_title)) {
      if (is.null(lay$annotations)) lay$annotations <- list()
      lay$annotations <- c(lay$annotations, list(
        list(
          text = this_title,
          x = mean(x_domains[[col_j]]),
          y = y_domains[[row_i]][2] + 0.01,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          font = list(size = 12),
          xanchor = "center"
        )
      ))
    }
  }
  
  # IMPORTANT: do.call() instead of !!!lay
  out <- do.call(plotly::layout, c(list(out), lay))
  out
}

multi_fig <- make_multiscene_3d(
  fig_list = list(fig_A, fig_B, fig_C, fig_D, fig_E, fig_F),
  nrows = 2, ncols = 3
)

outfile <- "FigureX_Conceptual_Geometry_AtoF_multiscene.html"
htmlwidgets::saveWidget(multi_fig, outfile, selfcontained = TRUE)

multi_fig

# -----------------------------
# 4b) Count convex-hull vertices per panel (A–F)
# -----------------------------
hull_vertex_count <- function(hull_obj) {
  stopifnot(!is.null(hull_obj$facets))
  length(unique(as.vector(hull_obj$facets)))
}

panel_vertex_counts <- c(
  A = hull_vertex_count(h_A),
  B = hull_vertex_count(h_B),
  C = hull_vertex_count(h_C),
  D = hull_vertex_count(h_D),
  E = hull_vertex_count(h_E),
  F = hull_vertex_count(h_F)
)

panel_n_points <- c(A=30, B=30, C=30, D=30, E=30, F=30)

panel_summary <- data.frame(
  Panel = names(panel_vertex_counts),
  n_points = as.integer(panel_n_points),
  n_hull_vertices = as.integer(panel_vertex_counts),
  stringsAsFactors = FALSE
)

print(panel_summary)

# -----------------------------
# 4c) Caption generator with vertex counts injected
# -----------------------------
make_figureX_caption <- function(panel_vertex_counts, n_points = 30L) {
  stopifnot(all(c("A","B","C","D","E","F") %in% names(panel_vertex_counts)))
  
  sprintf(
    "Figure X. Conceptual geometric representations of omic signature organization.
(A) Schematic illustration of an omic signature represented as a multidimensional informational entity projected into three-dimensional latent space. Individual points correspond to molecular components positioned according to their contributions across mechanistic, phenotypic, and contextual dimensions (n = %d points; convex-hull vertices = %d). The convex hull delineates the minimal geometric envelope enclosing the component cloud, providing an interpretable projection of higher-dimensional organization.
(B–C) Conceptual examples of distinct geometric fingerprints arising from different internal organization. Although the signatures contain comparable numbers of components (n = %d points each), differences in spatial arrangement yield compact, isotropic geometries (B; hull vertices = %d) or elongated, directional geometries (C; hull vertices = %d), reflecting coherent versus axis-dominated organization, respectively.
(D–F) Canonical geometric regimes illustrating organizational heterogeneity. Compact and isotropic hulls (D; hull vertices = %d) indicate coherent, single-axis behavior; elongated and anisotropic hulls (E; hull vertices = %d) reflect directional structure and partial divergence; irregular and faceted hulls (F; hull vertices = %d) indicate heterogeneous organization consistent with the coexistence of multiple semi-independent biological modules.
All geometries shown are generated from synthetic data and are intended solely to illustrate conceptual principles. These representations do not depict empirical results but serve as schematic guides for interpreting geometric organization in multidimensional omic signature space. (Source code and link to the corresponding .html figure are in Supporting Information).",
n_points, panel_vertex_counts["A"],
n_points, panel_vertex_counts["B"], panel_vertex_counts["C"],
panel_vertex_counts["D"], panel_vertex_counts["E"], panel_vertex_counts["F"]
  )
}

caption_FigureX <- make_figureX_caption(panel_vertex_counts, n_points = 30L)
cat(caption_FigureX)



###############################################################################
# SAVE FIGURE X PANELS (A–F) AS HTML + TIFF (600 dpi)
# Assumes fig_A ... fig_F already exist as plotly objects.
# FIX: removes magick::image_density() (not exported in some magick versions)
#      and embeds DPI via magick::image_write(..., density="600x600").
###############################################################################

suppressPackageStartupMessages({
  library(plotly)
  library(htmlwidgets)
})

# -----------------------------
# 0) Output folders
# -----------------------------
out_html <- "FigureX_Panels_HTML"
out_tiff <- "FigureX_Panels_TIFF_600dpi"
dir.create(out_html, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tiff, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Collect panels in a named list
# -----------------------------
fig_list <- list(
  A = fig_A,
  B = fig_B,
  C = fig_C,
  D = fig_D,
  E = fig_E,
  F = fig_F
)

# -----------------------------
# 2) HTML export (selfcontained)
# -----------------------------
for (nm in names(fig_list)) {
  f_html <- file.path(out_html, sprintf("FigureX_%s.html", nm))
  htmlwidgets::saveWidget(fig_list[[nm]], file = f_html, selfcontained = TRUE)
  message("Saved HTML: ", normalizePath(f_html))
}

# -----------------------------
# 3) TIFF 600 dpi export helpers
# -----------------------------

# Choose a pixel geometry that corresponds to your intended print size at 600 dpi.
# Example: 3.5 in × 3.0 in @ 600 dpi -> 2100 × 1800 px
tiff_width_in  <- 3.5
tiff_height_in <- 3.0
dpi_target     <- 600
px_w <- as.integer(round(tiff_width_in  * dpi_target))
px_h <- as.integer(round(tiff_height_in * dpi_target))

save_plotly_tiff_600dpi <- function(fig, out_file, width_px, height_px, dpi = 600) {
  stopifnot(inherits(fig, "plotly"))
  
  suppressPackageStartupMessages({
    library(htmlwidgets)
    library(webshot2)
    library(magick)
  })
  
  # Stable temp dir inside working directory (NOT /tmp)
  tmp_safe_dir <- file.path(getwd(), "_tmp_webshot")
  dir.create(tmp_safe_dir, showWarnings = FALSE, recursive = TRUE)
  
  stamp   <- paste0(Sys.getpid(), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample.int(1e6, 1))
  tmp_html <- file.path(tmp_safe_dir, paste0("plotly_", stamp, ".html"))
  tmp_png  <- file.path(tmp_safe_dir, paste0("plotly_", stamp, ".png"))
  
  # 1) Save widget to HTML (stable path)
  htmlwidgets::saveWidget(fig, file = tmp_html, selfcontained = TRUE)
  if (!file.exists(tmp_html)) stop("HTML was not created: ", tmp_html)
  
  # 2) Screenshot via Chromium (webshot2 expects a URL)
  webshot2::webshot(
    url     = paste0("file://", tmp_html),
    file    = tmp_png,
    vwidth  = width_px,
    vheight = height_px,
    zoom    = 1
  )
  if (!file.exists(tmp_png)) stop("PNG was not created: ", tmp_png)
  
  # 3) Convert to TIFF with embedded DPI
  img <- magick::image_read(tmp_png)
  img <- magick::image_resize(img, sprintf("%dx%d!", width_px, height_px))
  
  magick::image_write(
    img,
    path        = out_file,
    format      = "tiff",
    compression = "lzw",
    density     = paste0(dpi, "x", dpi)
  )
  
  if (!file.exists(out_file)) stop("TIFF was not created: ", out_file)
  invisible(out_file)
}

# -----------------------------
# 4) TIFF export (A–F)
# -----------------------------
for (nm in names(fig_list)) {
  f_tif <- file.path(out_tiff, sprintf("FigureX_%s_600dpi.tiff", nm))
  save_plotly_tiff_600dpi(
    fig      = fig_list[[nm]],
    out_file = f_tif,
    width_px = px_w,
    height_px= px_h,
    dpi      = dpi_target
  )
  message("Saved TIFF: ", normalizePath(f_tif))
}

###############################################################################
# DONE
###############################################################################


###
###
### Table 1
###
###
###############################################################################
# Distribution of circuitries by distance_implication
# Output: Excel (.xlsx) with counts + proportions (+ optional stratifications)
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(rio)
  library(readr)
})

# ---------------------------------------------------------------------------
# 0) Safety checks
# ---------------------------------------------------------------------------
stopifnot(exists("circuitries"))
stopifnot(is.data.frame(circuitries))

# Accept either "distance_implication" (correct) or the misspelled variant
dist_col <- dplyr::case_when(
  "distance_implication" %in% names(circuitries) ~ "distance_implication",
  "distace implication"  %in% names(circuitries) ~ "distace implication",
  TRUE ~ NA_character_
)

if (is.na(dist_col)) {
  stop("Column not found: expected 'distance_implication' (or 'distace implication').")
}

# ---------------------------------------------------------------------------
# 1) Global distribution (counts + proportions)
# ---------------------------------------------------------------------------
dist_global <- circuitries %>%
  mutate(distance_implication = .data[[dist_col]]) %>%
  mutate(distance_implication = ifelse(is.na(distance_implication) | distance_implication == "",
                                       "NA/blank", as.character(distance_implication))) %>%
  count(distance_implication, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(
    prop = n / sum(n),
    pct  = 100 * prop
  )
# ---------------------------------------------------------------------------
# 2) Optional: distribution by cancer-type prefix in Circuitries_id (e.g., "LGG")
#    (Keeps things tidy for downstream manuscript reporting)
# ---------------------------------------------------------------------------

dist_by_cancer <- circuitries %>%
  mutate(distance_implication = .data[[dist_col]]) %>%
  mutate(distance_implication = ifelse(is.na(distance_implication) | distance_implication == "",
                                       "NA/blank", as.character(distance_implication))) %>%
  mutate(
    cancer_type = sub("^([A-Za-z0-9]+).*", "\\1", as.character(.data$Circuitries_id))
  ) %>%
  filter(!is.na(cancer_type) & cancer_type != "") %>%
  count(cancer_type, distance_implication, name = "n") %>%
  group_by(cancer_type) %>%
  mutate(
    prop = n / sum(n),
    pct  = 100 * prop
  ) %>%
  ungroup() %>%
  arrange(cancer_type, desc(n))

# Wide-format pivot (useful as an Excel “matrix”)
dist_by_cancer_wide <- dist_by_cancer %>%
  dplyr::select(cancer_type, distance_implication, n) %>%
  tidyr::pivot_wider(
    names_from  = distance_implication,
    values_from = n,
    values_fill = 0
  )

# ---------------------------------------------------------------------------
# 3) Export to Excel (multiple sheets)
# ---------------------------------------------------------------------------
out_xlsx <- "Table 1. circuitries_distance_implication_distribution.xlsx"

rio::export(
  list(
    global_distribution          = dist_global,
    by_cancer_long_counts_props  = dist_by_cancer,
    by_cancer_wide_counts        = dist_by_cancer_wide
  ),
  file = out_xlsx
)

message("Saved Excel: ", normalizePath(out_xlsx))

####
####
#### Table 2
####
####
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)  # robust Excel writer
})

# ---------------------------------------------------------------------------
# 0) Safety checks
# ---------------------------------------------------------------------------
stopifnot(exists("circuitries"))
stopifnot(is.data.frame(circuitries))

# Accept either "vol_implication" (correct) or a misspelled variant
vol_col <- dplyr::case_when(
  "vol_implication" %in% names(circuitries) ~ "vol_implication",
  "vol implication"  %in% names(circuitries) ~ "vol implication",
  TRUE ~ NA_character_
)
if (is.na(vol_col)) {
  stop("Column not found: expected 'vol_implication' (or 'vol implication').")
}

has_circuitry_id <- "Circuitries_id" %in% names(circuitries)

# ---------------------------------------------------------------------------
# 1) Global distribution: counts + proportions + percentages
# ---------------------------------------------------------------------------
vol_global <- circuitries %>%
  mutate(vol_implication = .data[[vol_col]]) %>%
  mutate(vol_implication = ifelse(is.na(vol_implication) | vol_implication == "",
                                  "NA/blank", as.character(vol_implication))) %>%
  count(vol_implication, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(
    prop = n / sum(n),
    pct  = 100 * prop
  )

# Nice formatting for reporting (optional)
vol_global_report <- vol_global %>%
  mutate(
    prop = round(prop, 6),
    pct  = round(pct, 1)
  ) %>%
  rename(`Volume-implication regime` = vol_implication,
         `Proportion` = prop,
         `Percentage (%)` = pct)

# ---------------------------------------------------------------------------
# 2) Optional: distribution by cancer-type prefix (from Circuitries_id)
# ---------------------------------------------------------------------------
vol_by_cancer_report <- NULL

if (has_circuitry_id) {
  vol_by_cancer <- circuitries %>%
    mutate(vol_implication = .data[[vol_col]]) %>%
    mutate(vol_implication = ifelse(is.na(vol_implication) | vol_implication == "",
                                    "NA/blank", as.character(vol_implication))) %>%
    mutate(
      cancer_type = sub("^([A-Za-z0-9]+).*", "\\1", as.character(Circuitries_id))
    ) %>%
    count(cancer_type, vol_implication, name = "n") %>%
    group_by(cancer_type) %>%
    mutate(
      prop = n / sum(n),
      pct  = 100 * prop
    ) %>%
    ungroup() %>%
    arrange(cancer_type, desc(n))
  
  vol_by_cancer_report <- vol_by_cancer %>%
    mutate(
      prop = round(prop, 6),
      pct  = round(pct, 1)
    ) %>%
    rename(`Cancer type` = cancer_type,
           `Volume-implication regime` = vol_implication,
           `Proportion` = prop,
           `Percentage (%)` = pct)
}

# ---------------------------------------------------------------------------
# 3) Export to Excel (single workbook with two sheets)
# ---------------------------------------------------------------------------
out_xlsx <- "Table 2. circuitries_vol_implication_distribution.xlsx"

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "vol_implication_global")
openxlsx::writeData(wb, "vol_implication_global", vol_global_report)

if (!is.null(vol_by_cancer_report)) {
  openxlsx::addWorksheet(wb, "vol_implication_by_cancer")
  openxlsx::writeData(wb, "vol_implication_by_cancer", vol_by_cancer_report)
}

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Saved Excel workbook: ", normalizePath(out_xlsx))

# Return objects for immediate inspection
vol_global_report

###############################################################################
# Volume-implication regime distribution
# Manuscript-consistent, Excel-ready
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

# ---------------------------------------------------------------------------
# 0) Safety checks
# ---------------------------------------------------------------------------
stopifnot(exists("circuitries"))
stopifnot(is.data.frame(circuitries))
stopifnot("vol_implication" %in% names(circuitries))

# ---------------------------------------------------------------------------
# 1) Manuscript-consistent label mapping
#    (DO NOT change wording – aligned with Methods/Results)
# ---------------------------------------------------------------------------
vol_labels <- c(
  "low_dimensional_flat ; symmetric_sig_int" =
    "Low-dimensional, symmetric geometry",
  "low_dimensional_flat ; moderately_asymmetric_sig_int" =
    "Low-dimensional, moderately asymmetric geometry",
  "low_dimensional_flat ; strongly_asymmetric_sig_int" =
    "Low-dimensional, strongly asymmetric geometry",
  
  "intermediate_complexity ; symmetric_sig_int" =
    "Intermediate-complexity, symmetric geometry",
  "intermediate_complexity ; moderately_asymmetric_sig_int" =
    "Intermediate-complexity, moderately asymmetric geometry",
  "intermediate_complexity ; strongly_asymmetric_sig_int" =
    "Intermediate-complexity, strongly asymmetric geometry",
  
  "high_complexity_multidimensional ; symmetric_sig_int" =
    "High-complexity multidimensional, symmetric geometry",
  "high_complexity_multidimensional ; moderately_asymmetric_sig_int" =
    "High-complexity multidimensional, moderately asymmetric geometry",
  "high_complexity_multidimensional ; strongly_asymmetric_sig_int" =
    "High-complexity multidimensional, strongly asymmetric geometry"
)

# ---------------------------------------------------------------------------
# 2) Build distribution table
# ---------------------------------------------------------------------------
vol_distribution <- circuitries %>%
  mutate(
    Volume_implication_internal = vol_implication,
    Geometric_volume_implication_regime =
      unname(vol_labels[Volume_implication_internal])
  ) %>%
  count(Geometric_volume_implication_regime, name = "Count (n)") %>%
  arrange(desc(`Count (n)`)) %>%
  mutate(
    Proportion      = `Count (n)` / sum(`Count (n)`),
    `Percentage (%)` = round(100 * Proportion, 1)
  )

# ---------------------------------------------------------------------------
# 3) Inspect in R (optional but recommended)
# ---------------------------------------------------------------------------
print(vol_distribution)

# ---------------------------------------------------------------------------
# 4) Export to Excel (journal-safe)
# ---------------------------------------------------------------------------
out_xlsx <- "Table_Volume_Implication_Distribution.xlsx"

wb <- createWorkbook()
addWorksheet(wb, "volume_implication_distribution")
writeData(wb, "volume_implication_distribution", vol_distribution)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Excel table saved at: ", normalizePath(out_xlsx))

###############################################################################
# END
###############################################################################

###############################################################################
# Figure 4 exemplar picker + builder + saver (A–D)
# - Keeps A/B fixed (user-specified)
# - Auto-selects C/D: distinct, non-duplicate (order-invariant), non-collapsed
# - Prints an explicit report: which circuitries were used + where files saved
# - Saves 4 self-contained HTML plotly panels into out_dir
###############################################################################
make_figure4_exemplars <- function(circuitries, tensor, embedding,
                                   A_id = "LGG-5422 / LGG-6519",
                                   B_id = "LGG-6708 / LGG-5904",
                                   out_dir = "Figure_4_panels_html",
                                   seed = 1,
                                   # collapse filters (used only if those columns exist)
                                   min_bary_dist = 0.15,
                                   min_hull_vol  = 1e-6,
                                   # plot save size
                                   width = 1600,
                                   height = 1200) {
  
  # ---- basic checks ----------------------------------------------------------
  stopifnot(is.data.frame(circuitries))
  stopifnot("Circuitries_id" %in% names(circuitries))
  
  # ---- packages --------------------------------------------------------------
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("stringr", quietly = TRUE))
  stopifnot(requireNamespace("htmlwidgets", quietly = TRUE))
  
  # ---- helpers ---------------------------------------------------------------
  # Normalize "X / Y" so "X / Y" == "Y / X" (order-invariant)
  normalize_id <- function(x) {
    x <- as.character(x)
    parts <- stringr::str_split(x, "\\s*/\\s*", simplify = TRUE)
    if (ncol(parts) < 2) return(x)
    a <- parts[, 1]; b <- parts[, 2]
    ifelse(a <= b, paste0(a, " / ", b), paste0(b, " / ", a))
  }
  
  # Safer filename tag
  sanitize_for_filename <- function(x) {
    x |>
      stringr::str_replace_all("\\s*/\\s*", "__") |>
      stringr::str_replace_all("[^A-Za-z0-9_\\-\\.]+", "_")
  }
  
  # ---- normalize + de-duplicate ---------------------------------------------
  # Keep one canonical original Circuitries_id per normalized key
  tbl <- circuitries |>
    dplyr::mutate(
      Circuitries_id = as.character(.data$Circuitries_id),
      circuitry_key  = normalize_id(.data$Circuitries_id)
    ) |>
    dplyr::group_by(.data$circuitry_key) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  
  A_key <- normalize_id(A_id)
  B_key <- normalize_id(B_id)
  
  if (!(A_key %in% tbl$circuitry_key)) stop("A_id not found in circuitries: ", A_id)
  if (!(B_key %in% tbl$circuitry_key)) stop("B_id not found in circuitries: ", B_id)
  
  # Candidate pool excluding A/B
  cand <- tbl |>
    dplyr::filter(!(.data$circuitry_key %in% c(A_key, B_key)))
  
  # Optional filters (only applied if the columns exist)
  if ("barycenter_distance" %in% names(cand)) {
    cand <- cand |>
      dplyr::filter(.data$barycenter_distance >= min_bary_dist)
  }
  if (all(c("sig_hull_vol", "int_hull_vol") %in% names(cand))) {
    cand <- cand |>
      dplyr::filter(.data$sig_hull_vol >= min_hull_vol,
                    .data$int_hull_vol >= min_hull_vol)
  }
  
  if (nrow(cand) < 2) {
    stop(
      "Not enough candidates for C/D after filters.\n",
      "Try lowering min_bary_dist/min_hull_vol, or verify circuitries table columns."
    )
  }
  
  set.seed(seed)
  pick_keys <- sample(cand$circuitry_key, size = 2, replace = FALSE)
  
  exemplar_keys <- c(A = A_key, B = B_key, C = pick_keys[[1]], D = pick_keys[[2]])
  
  # Map keys -> canonical original IDs (so builder receives an ID it recognizes)
  key_to_id <- setNames(tbl$Circuitries_id, tbl$circuitry_key)
  
  # ---- PATCH: keep *exact* user-provided orientation for A/B if present ------
  # This prevents "LGG-6708 / LGG-5904" being replaced by "LGG-5904 / LGG-6708"
  # when both exist but tbl kept the first occurrence.
  if (A_id %in% circuitries$Circuitries_id) key_to_id[[A_key]] <- A_id
  if (B_id %in% circuitries$Circuitries_id) key_to_id[[B_key]] <- B_id
  
  exemplar_ids <- unname(key_to_id[exemplar_keys])
  names(exemplar_ids) <- names(exemplar_keys)
  
  # ---- output directory ------------------------------------------------------
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_dir_abs <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  
  # ---- builder + saver -------------------------------------------------------
  build_and_save_panel <- function(panel_label, Circuitries_id) {
    
    # IMPORTANT: build_circuitry_polytope expects 'circuitry_id'
    poly <- build_circuitry_polytope(
      tensor = tensor,
      embedding = embedding,
      circuitry_id = Circuitries_id
    )
    
    # Your plotting function should return a plotly/htmlwidget
    p <- plot_circuitry_polytope_custom(
      poly_obj = poly,
      title = paste0("Panel ", panel_label, " — ", Circuitries_id),
      outfile_html = NULL
    )
    
    out_file <- file.path(
      out_dir_abs,
      paste0("Figure4_", panel_label, "_",
             sanitize_for_filename(Circuitries_id),
             "_dual_polytope_18D.html")
    )
    
    htmlwidgets::saveWidget(p, file = out_file, selfcontained = TRUE)
    
    message("Saved Panel ", panel_label, ": ", out_file)
    out_file
  }
  
  saved_files <- vapply(names(exemplar_ids), function(lbl) {
    build_and_save_panel(lbl, exemplar_ids[[lbl]])
  }, character(1))
  
  # ---- explicit reporting ----------------------------------------------------
  report <- data.frame(
    panel = names(exemplar_ids),
    Circuitries_id = unname(exemplar_ids),
    circuitry_key = unname(exemplar_keys),
    html_file = unname(saved_files),
    stringsAsFactors = FALSE
  )
  
  cat("\n================ Figure 4 exemplar selection ================\n")
  print(report, row.names = FALSE)
  cat("\nSaved HTML directory:\n  ", out_dir_abs, "\n", sep = "")
  
  invisible(report)
}

# Create an object in the environment with the selected Figure 4 circuitries
figure4_exemplars_2 <- make_figure4_exemplars(
  circuitries = circuitries,
  tensor = tensor,
  embedding = embedding,
  A_id = "LGG-3207 / LGG-2427",
  B_id = "BRCA-8458 / BRCA-2923",
  out_dir = "figures_html",
  seed = 7
)

# Inspect it
figure4_exemplars_2

###############################################################################
# RUN (example)
###############################################################################
# Assumes you already have:
#   - circuitries (data.frame with column 'Circuitries_id', and optionally
#       'barycenter_distance', 'sig_hull_vol', 'int_hull_vol')
#   - tensor, embedding in memory
#   - functions: build_circuitry_polytope(), plot_circuitry_polytope_custom()

report <- make_figure4_exemplars(
  circuitries = circuitries,
  tensor = tensor,
  embedding = embedding,
  A_id = "LGG-5422 / LGG-6519",
  B_id = "LGG-6708 / LGG-5904",
  out_dir = "Figure_4_panels_html",
  seed = 7,
  min_bary_dist = 0.3,
  min_hull_vol  = 0.01
)

report

figure4_exemplars_2 <- report

###############################################################################
# EXTRACT RELEVANT VARIABLES FOR THE CIRCUITRY: LUSC-6681 / LUSC-8604
###############################################################################

# 1. Target circuitry identifier
target_id <- "LUSC-6681 / LUSC-8604"

# 2. Variables to extract (same structure as cir_relevant_variables.tsv)
vars_to_extract <- c(
  "Circuitries_id",
  "Nomenclature_sig",
  "Signatures",
  "Nomenclature_int",
  "Interaction",
  "barycenter_distance",
  "sig_hull_vol",
  "int_hull_vol",
  "vol_ratio",
  "distance_implication",
  "vol_implication",
  "Omic_layer_sig",
  "Phenotypic_layer_sig",
  "Omic_layer_int",
  "Phenotypic_layer_int",
  "Correlation_rho_sig",
  "Correlation_rho_int",
  "Tumor_vs_normal_sig",
  "Tumor_vs_normal_int",
  "Combined_outcome_HRC_sig",
  "Combined_outcome_HRC_int",
  "Microenvironment_classification_sig",
  "Microenvironment_classification_int",
  "Immune_classification_sig",
  "Immune_classification_int",
  "Phenotypic_concordance",
  "Survival_concordance_aggregated",
  "Immune_concordance",
  "Final_concordance_summary",
  "CTAB",
  "Metabolism",
  "Pathways"
)
# 3. Extract, convert all values to character (required for tidy pivot), pivot to two-column format
cir_6681_8604_relevant_variables <- circuitries %>%
  dplyr::filter(Circuitries_id == target_id) %>%
  dplyr::select(dplyr::all_of(vars_to_extract)) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%  # ensure homogeneous type
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# 4. Print result
cir_6681_8604_relevant_variables

# 5. (Optional) Save to TSV
readr::write_tsv(cir_6681_8604_relevant_variables, "LUSC_6681_LUSC_8604_cir_relevant_variables.tsv")
###############################################################################

###############################################################################
# EXTRACT RELEVANT VARIABLES FOR THE CIRCUITRY: LUAD-10062 / LUAD-11695
###############################################################################

# 1. Target circuitry identifier
target_id <- "LUAD-10062 / LUAD-11695"

# 2. Variables to extract (same structure as cir_relevant_variables.tsv)
vars_to_extract <- c(
  "Circuitries_id",
  "Nomenclature_sig",
  "Signatures",
  "Nomenclature_int",
  "Interaction",
  "barycenter_distance",
  "sig_hull_vol",
  "int_hull_vol",
  "vol_ratio",
  "distance_implication",
  "vol_implication",
  "Omic_layer_sig",
  "Phenotypic_layer_sig",
  "Omic_layer_int",
  "Phenotypic_layer_int",
  "Correlation_rho_sig",
  "Correlation_rho_int",
  "Tumor_vs_normal_sig",
  "Tumor_vs_normal_int",
  "Combined_outcome_HRC_sig",
  "Combined_outcome_HRC_int",
  "Microenvironment_classification_sig",
  "Microenvironment_classification_int",
  "Immune_classification_sig",
  "Immune_classification_int",
  "Phenotypic_concordance",
  "Survival_concordance_aggregated",
  "Immune_concordance",
  "Final_concordance_summary",
  "CTAB",
  "Metabolism",
  "Pathways"
)
# 3. Extract, convert all values to character (required for tidy pivot), pivot to two-column format
cir_10062_11695_relevant_variables <- circuitries %>%
  dplyr::filter(Circuitries_id == target_id) %>%
  dplyr::select(dplyr::all_of(vars_to_extract)) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%  # ensure homogeneous type
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# 4. Print result
cir_10062_11695_relevant_variables

# 5. (Optional) Save to TSV
readr::write_tsv(cir_10062_11695_relevant_variables, "LUAD_10062_LUAD_11695_cir_relevant_variables.tsv")
###############################################################################


###############################################################################
# EXTRACT RELEVANT VARIABLES FOR THE CIRCUITRY: LGG-5422 / LGG-6519
###############################################################################

# 1. Target circuitry identifier
target_id <- "LGG-5422 / LGG-6519"

# 2. Variables to extract (same structure as cir_relevant_variables.tsv)
vars_to_extract <- c(
  "Circuitries_id",
  "Nomenclature_sig",
  "Signatures",
  "Nomenclature_int",
  "Interaction",
  "barycenter_distance",
  "sig_hull_vol",
  "int_hull_vol",
  "vol_ratio",
  "distance_implication",
  "vol_implication",
  "Omic_layer_sig",
  "Phenotypic_layer_sig",
  "Omic_layer_int",
  "Phenotypic_layer_int",
  "Correlation_rho_sig",
  "Correlation_rho_int",
  "Tumor_vs_normal_sig",
  "Tumor_vs_normal_int",
  "Combined_outcome_HRC_sig",
  "Combined_outcome_HRC_int",
  "Microenvironment_classification_sig",
  "Microenvironment_classification_int",
  "Immune_classification_sig",
  "Immune_classification_int",
  "Phenotypic_concordance",
  "Survival_concordance_aggregated",
  "Immune_concordance",
  "Final_concordance_summary",
  "CTAB",
  "Metabolism",
  "Pathways"
)

# 3. Extract, convert all values to character (required for tidy pivot), pivot to two-column format
cir_5422_6519_relevant_variables <- circuitries %>%
  dplyr::filter(Circuitries_id == target_id) %>%
  dplyr::select(dplyr::all_of(vars_to_extract)) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%  # ensure homogeneous type
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# 4. Print result
cir_5422_6519_relevant_variables

# 5. (Optional) Save to TSV
readr::write_tsv(cir_5422_6519_relevant_variables, "LGG_5422_LGG_6519_cir_relevant_variables.tsv")
###############################################################################

###############################################################################
# Filter 'circuitries' for a fixed set of Circuitries_id (order-invariant)
###############################################################################

stopifnot(exists("circuitries"))
stopifnot(is.data.frame(circuitries))
stopifnot("Circuitries_id" %in% names(circuitries))

# ---- target circuitries (as provided) ---------------------------------------
target_ids <- c(
  "LGG-5422 / LGG-6519",
  "LGG-6708 / LGG-5904",
  "LUSC-6681 / LUSC-8604",
  "LUAD-10062 / LUAD-11695"
)

# ---- normalize "X / Y" so "X / Y" == "Y / X" --------------------------------
normalize_id <- function(x) {
  x <- as.character(x)
  parts <- strsplit(x, "\\s*/\\s*")
  vapply(parts, function(p) {
    if (length(p) < 2) return(x)
    a <- p[1]; b <- p[2]
    if (a <= b) paste0(a, " / ", b) else paste0(b, " / ", a)
  }, character(1))
}

# ---- filter -----------------------------------------------------------------
target_keys <- normalize_id(target_ids)

circuitries_F4 <- circuitries[
  normalize_id(circuitries$Circuitries_id) %in% target_keys,
  ,
  drop = FALSE
]

# Optional: inspect
circuitries_F4$Circuitries_id
nrow(circuitries_F4)

stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("stringr", quietly = TRUE))
stopifnot("Circuitries_id" %in% names(circuitries))

###############################################################################
# AUDIT + NORMALIZE + DEDUPLICATE Circuitries_id
# Goal:
#  (1) Strictly audit whether Circuitries_id follows EXACT canonical format:
#      "CANCER-NNNN / CANCER-NNNN"  (ONE space, slash, ONE space)
#  (2) Detect whitespace/pathology variants (leading/trailing, multi-spaces, etc.)
#  (3) Detect redundant circuitries due to swapped order (A / B vs B / A)
#  (4) Provide an EXACT-MATCH filter that returns ONLY the 4 requested rows
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

stopifnot(exists("circuitries"))
stopifnot(is.data.frame(circuitries))
stopifnot("Circuitries_id" %in% names(circuitries))

# -----------------------------------------------------------------------------
# 1) Define canonical strict pattern
#    EXACTLY: token-digits␠/␠token-digits
# -----------------------------------------------------------------------------
canon_pat <- "^[A-Za-z0-9]+-[0-9]+ / [A-Za-z0-9]+-[0-9]+$"

# -----------------------------------------------------------------------------
# 2) Core helpers
# -----------------------------------------------------------------------------
trim_ws <- function(x) str_trim(as.character(x), side = "both")

# Collapse any whitespace into single spaces, then enforce " / " around slash
# This is a *normalizer*; use only if you want to repair IDs, not just audit.
normalize_circuitry_id_strict <- function(x) {
  x <- trim_ws(x)
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_replace_all(x, "\\s*/\\s*", " / ")
  x
}

# Order-invariant key so "A / B" and "B / A" map to the same identity
# (Assumes normalized delimiter " / " is present)
make_order_invariant_key <- function(x) {
  x <- normalize_circuitry_id_strict(x)
  parts <- str_split_fixed(x, " / ", 2)
  # If split fails, return as-is (will be flagged elsewhere)
  if (ncol(parts) < 2) return(x)
  a <- parts[, 1]
  b <- parts[, 2]
  ifelse(a <= b, paste0(a, " / ", b), paste0(b, " / ", a))
}

# -----------------------------------------------------------------------------
# 3) AUDIT TABLE: identify format violations and whitespace anomalies
# -----------------------------------------------------------------------------
audit_tbl <- circuitries %>%
  transmute(
    Circuitries_id_raw = as.character(.data$Circuitries_id),
    Circuitries_id_trim = trim_ws(.data$Circuitries_id),
    Circuitries_id_norm = normalize_circuitry_id_strict(.data$Circuitries_id),
    
    # STRICT format checks
    is_canonical_raw  = str_detect(.data$Circuitries_id_raw,  canon_pat),
    is_canonical_trim = str_detect(.data$Circuitries_id_trim, canon_pat),
    is_canonical_norm = str_detect(.data$Circuitries_id_norm, canon_pat),
    
    # Common pathology flags
    has_leading_or_trailing_ws = (.data$Circuitries_id_raw != .data$Circuitries_id_trim),
    has_any_tab_or_newline     = str_detect(.data$Circuitries_id_raw, "[\\t\\r\\n]"),
    has_multiple_spaces        = str_detect(.data$Circuitries_id_raw, " {2,}"),
    has_no_spaces_around_slash = str_detect(.data$Circuitries_id_raw, "/") &
      !str_detect(.data$Circuitries_id_raw, " / "),
    has_weird_slash_spacing    = str_detect(.data$Circuitries_id_raw, "\\s*/\\s*") &
      !str_detect(.data$Circuitries_id_raw, " / ")
  )

# Summary counts
audit_summary <- audit_tbl %>%
  summarise(
    n_total = n(),
    n_canonical_raw  = sum(is_canonical_raw,  na.rm = TRUE),
    n_canonical_trim = sum(is_canonical_trim, na.rm = TRUE),
    n_canonical_norm = sum(is_canonical_norm, na.rm = TRUE),
    
    n_leading_trailing_ws = sum(has_leading_or_trailing_ws, na.rm = TRUE),
    n_tabs_newlines       = sum(has_any_tab_or_newline,     na.rm = TRUE),
    n_multi_spaces        = sum(has_multiple_spaces,        na.rm = TRUE),
    n_no_spaces_around_sl = sum(has_no_spaces_around_slash, na.rm = TRUE),
    n_weird_slash_spacing = sum(has_weird_slash_spacing,    na.rm = TRUE)
  )

cat("\n================ Circuitries_id AUDIT SUMMARY ================\n")
print(audit_summary)

# Show a small set of non-canonical raw examples (if any)
noncanon_examples <- audit_tbl %>%
  filter(!is_canonical_raw) %>%
  distinct(Circuitries_id_raw) %>%
  slice_head(n = 20)

if (nrow(noncanon_examples) > 0) {
  cat("\n--- Examples of NON-canonical Circuitries_id (raw) ---\n")
  print(noncanon_examples, row.names = FALSE)
} else {
  cat("\nAll Circuitries_id values are canonical under the STRICT raw pattern.\n")
}

# -----------------------------------------------------------------------------
# 4) REDUNDANCY CHECK: swapped-order duplicates
# -----------------------------------------------------------------------------
dup_tbl <- circuitries %>%
  mutate(
    Circuitries_id_norm = normalize_circuitry_id_strict(.data$Circuitries_id),
    circuitry_key = make_order_invariant_key(.data$Circuitries_id)
  ) %>%
  count(circuitry_key, name = "n_rows") %>%
  filter(n_rows > 1) %>%
  arrange(desc(n_rows))

cat("\n================ Swapped-order DUPLICATES (key-level) ================\n")
if (nrow(dup_tbl) == 0) {
  cat("No swapped-order duplicates detected (order-invariant keys are unique).\n")
} else {
  print(dup_tbl, row.names = FALSE)
}

# OPTIONAL: inspect the first duplicated key’s actual raw rows
if (nrow(dup_tbl) > 0) {
  first_key <- dup_tbl$circuitry_key[1]
  cat("\n--- Raw rows for first duplicated key ---\n")
  print(
    circuitries %>%
      mutate(
        Circuitries_id_norm = normalize_circuitry_id_strict(.data$Circuitries_id),
        circuitry_key = make_order_invariant_key(.data$Circuitries_id)
      ) %>%
      filter(circuitry_key == first_key) %>%
      select(Circuitries_id, Circuitries_id_norm, circuitry_key) %>%
      distinct(),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 5) EXACT-MATCH FILTER (STRICT): returns ONLY rows whose *raw* string matches
#    exactly one of your provided IDs (including spaces).
#    This will NOT match swapped order.
# -----------------------------------------------------------------------------
target_ids_exact <- c(
  "LGG-5422 / LGG-6519",
  "LGG-6708 / LGG-5904",
  "LUSC-6681 / LUSC-8604",
  "LUAD-10062 / LUAD-11695"
)

circuitries_exact4 <- circuitries %>%
  filter(.data$Circuitries_id %in% target_ids_exact)

cat("\n================ EXACT-MATCH FILTER RESULT ================\n")
cat("Rows returned: ", nrow(circuitries_exact4), "\n", sep = "")
print(circuitries_exact4 %>% select(Circuitries_id) %>% distinct(), row.names = FALSE)

# -----------------------------------------------------------------------------
# 6) ORDER-INVARIANT FILTER (OPTIONAL): matches either orientation by key.
#    Use this if you want to retrieve rows regardless of "A / B" vs "B / A".
# -----------------------------------------------------------------------------
target_keys <- make_order_invariant_key(target_ids_exact)

circuitries_by_key <- circuitries %>%
  mutate(circuitry_key = make_order_invariant_key(.data$Circuitries_id)) %>%
  filter(circuitry_key %in% target_keys)

cat("\n================ ORDER-INVARIANT FILTER RESULT (optional) ================\n")
cat("Rows returned: ", nrow(circuitries_by_key), "\n", sep = "")
print(
  circuitries_by_key %>%
    select(Circuitries_id, circuitry_key) %>%
    distinct() %>%
    arrange(circuitry_key, Circuitries_id),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5) Figure 4 IDs EXACT-MATCH FILTER (STRICT): returns ONLY rows whose *raw* string matches
#    exactly one of your provided IDs (including spaces).
#    This will NOT match swapped order.
# -----------------------------------------------------------------------------
target_ids_exact <- c(
  "LGG-5422 / LGG-6519",
  "LGG-6708 / LGG-5904",
  "LUSC-6681 / LUSC-8604",
  "LUAD-10062 / LUAD-11695"
)

circuitries_exact4 <- circuitries %>%
  filter(.data$Circuitries_id %in% target_ids_exact)

cat("\n================ EXACT-MATCH FILTER RESULT ================\n")
cat("Rows returned: ", nrow(circuitries_exact4), "\n", sep = "")
print(circuitries_exact4 %>% select(Circuitries_id) %>% distinct(), row.names = FALSE)

# -----------------------------------------------------------------------------
# OPTIONAL BUT RECOMMENDED: expanded raw rows for ALL duplicated keys
#     - This prevents the “it became smaller” confusion: this table keeps ALL rows
# -----------------------------------------------------------------------------
dup_rows <- circuitries_original %>%
  dplyr::mutate(
    Circuitries_id_norm = normalize_circuitry_id_strict(.data$Circuitries_id),
    circuitry_key       = make_order_invariant_key(.data$Circuitries_id)
  ) %>%
  dplyr::semi_join(dup_tbl %>% dplyr::select(.data$circuitry_key), by = "circuitry_key") %>%
  dplyr::arrange(.data$circuitry_key, .data$Circuitries_id_norm)

cat("\n--- Raw rows belonging to duplicated keys (dup_rows) ---\n")
cat("nrow(dup_rows) = ", nrow(dup_rows), "\n", sep = "")

# ------------------------------------------------------------------------------
# Snippet to generate and annotate swapped regulatory circuitries
# ------------------------------------------------------------------------------

library(dplyr)

circuitries_annotated <- circuitries %>%
  mutate(
    # (1) Observed circuitry identifier, strictly normalized
    #     This preserves the original ordering as reported in the dataset
    circuitry_id_normalized = normalize_circuitry_id_strict(.data$Circuitries_id),
    
    # (2) Canonical biological identity (order-invariant)
    #     This collapses equivalent circuitries regardless of signature–interaction order
    circuitry_id_canonical  = make_order_invariant_key(.data$Circuitries_id),
    
    # (3) Observed orientation relative to the canonical representation
    #     Indicates whether the circuitry appears in canonical or swapped order
    circuitry_orientation = if_else(
      circuitry_id_normalized == circuitry_id_canonical,
      "canonical_order",
      "swapped_order"
    )
  ) %>%
  group_by(.data$circuitry_id_canonical) %>%
  mutate(
    # (4) Final flag reported in the manuscript
    #     Marks circuitries that appear more than once with opposite orientations
    swapped_circuitry = if_else(n() > 1, "Yes", "No")
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Export annotated circuitry tables
# ------------------------------------------------------------------------------

rio::export(circuitries_annotated, "Regulatory_circuitries_02.tsv")
rio::export(circuitries_annotated, "Dataset_S5_amended.tsv")
rio::export(circuitries_annotated, "Dataset_S5_amended.rds")

###############################################################################
# Figure 4: TRUE 2×2 composite WITHOUT merging plotly objects
# - Uses the 4 already-generated standalone HTML panels (A–D)
# - Creates ONE master HTML with a CSS grid (2×2) of iframes
# - Absolutely no trace/scene overlap possible (each panel isolated)
###############################################################################

stopifnot(requireNamespace("htmltools", quietly = TRUE))

# Must exist from your earlier step:
stopifnot(exists("report"))
stopifnot(all(c("panel", "html_file") %in% names(report)))

# ---- ensure A–D order -------------------------------------------------------
panel_order <- c("A","B","C","D")
tbl <- report
tbl$panel <- as.character(tbl$panel)
tbl <- tbl[order(match(tbl$panel, panel_order)), ]

# ---- debug: Check file paths -------------------------------------------------
cat("Checking HTML file paths:\n")
cat(tbl$html_file, sep = "\n")

# ---- validate files ---------------------------------------------------------
# Ensure all HTML files exist
tbl$html_file <- normalizePath(tbl$html_file, winslash = "/", mustWork = FALSE)

# Debugging check for existence of the files
file_exists <- file.exists(tbl$html_file)
if (any(!file_exists)) {
  cat("Some files do not exist:\n")
  cat(tbl$html_file[!file_exists], sep = "\n")
  stop("Please check that all files exist and are correct.")
}

# ---- put composite in SAME folder as the panels -----------------------------
out_dir <- normalizePath(dirname(tbl$html_file[1]), winslash = "/", mustWork = TRUE)

# If panels are not all in the same folder, copy them into out_dir
# (So the iframe src can be just the basename)
if (length(unique(dirname(normalizePath(tbl$html_file, winslash="/")))) != 1) {
  message("Panels are in different folders; copying all panels into: ", out_dir)
  for (i in seq_len(nrow(tbl))) {
    src <- normalizePath(tbl$html_file[i], winslash = "/", mustWork = TRUE)
    dst <- file.path(out_dir, basename(src))
    if (!file.exists(dst)) file.copy(src, dst, overwrite = TRUE)
    tbl$html_file[i] <- dst
  }
}

# ---- build master HTML ------------------------------------------------------
# Tweak iframe height if you want taller panels:
iframe_height_px <- 900

make_iframe <- function(src, title) {
  htmltools::tags$iframe(
    src   = basename(src),   # relative to master html in same folder
    style = paste0(
      "width:100%; height:", iframe_height_px, "px; ",
      "border:0; border-radius:10px; ",
      "box-shadow:0 2px 12px rgba(0,0,0,0.12);"
    ),
    `aria-label` = title
  )
}

page <- htmltools::tags$html(
  htmltools::tags$head(
    htmltools::tags$meta(charset = "utf-8"),
    htmltools::tags$title("Figure 4 — Composite 2×2"),
    htmltools::tags$style(htmltools::HTML("
      body { margin: 0; padding: 18px; font-family: sans-serif; background: #ffffff; }
      .grid {
        display: grid;
        grid-template-columns: 1fr 1fr;
        grid-auto-rows: auto;
        gap: 18px;
        align-items: stretch;
      }
      .cell { width: 100%; }
      @media (max-width: 1200px) {
        .grid { grid-template-columns: 1fr; }
      }
    "))
  ),
  htmltools::tags$body(
    htmltools::tags$div(
      class = "grid",
      htmltools::tags$div(class="cell", make_iframe(tbl$html_file[tbl$panel=="A"], "Panel A")),
      htmltools::tags$div(class="cell", make_iframe(tbl$html_file[tbl$panel=="B"], "Panel B")),
      htmltools::tags$div(class="cell", make_iframe(tbl$html_file[tbl$panel=="C"], "Panel C")),
      htmltools::tags$div(class="cell", make_iframe(tbl$html_file[tbl$panel=="D"], "Panel D"))
    )
  )
)

out_file <- file.path(out_dir, "Figure4_COMPOSITE_2x2_iframe.html")
htmltools::save_html(page, file = out_file)

message("Saved TRUE 2×2 composite (no overlap possible): ", out_file)

# --- robust open (Windows) ---------------------------------------------------
out_file_abs <- normalizePath(out_file, winslash = "\\", mustWork = TRUE)
message("Opening: ", out_file_abs)

if (.Platform$OS.type == "windows") {
  shell.exec(out_file_abs)   # most reliable on Windows
} else {
  browseURL(out_file_abs)
}


#####
#####
##### Figure 4 Final
#####
#####
###############################################################################
# Figure 4 Final — Build 4 plotly panels (A–D) and assemble TRUE 2×2 multiscene HTML
# + Two-tier right-side legend WITHOUT overlap:
#   (1) TOP: Latent dimensions ((axes) 18D) annotation block
#   (2) BELOW: single Plotly trace legend (vertices/barycenters/hulls) shown once
###############################################################################

stopifnot(requireNamespace("plotly", quietly = TRUE))
stopifnot(requireNamespace("htmlwidgets", quietly = TRUE))

# ---------------------------------------------------------------------------
# 0) Inputs assumed to exist in your environment:
#   - report (data.frame) with columns: panel (A–D), Circuitries_id
#   - tensor, embedding
#   - functions: build_circuitry_polytope(), plot_circuitry_polytope_custom()
# ---------------------------------------------------------------------------

stopifnot(exists("report"))
stopifnot(all(c("panel","Circuitries_id") %in% names(report)))
stopifnot(exists("tensor"))
stopifnot(exists("embedding"))
stopifnot(exists("build_circuitry_polytope"))
stopifnot(exists("plot_circuitry_polytope_custom"))

panel_order <- c("A","B","C","D")
tbl <- report
tbl$panel <- as.character(tbl$panel)
tbl <- tbl[order(match(tbl$panel, panel_order)), ]
stopifnot(nrow(tbl) == 4)
stopifnot(all(tbl$panel %in% panel_order))

# ---------------------------------------------------------------------------
# 1) Build the FOUR plotly objects in memory
# ---------------------------------------------------------------------------

build_panel_plotly <- function(panel_label, circuitry_id, tensor, embedding) {
  
  poly <- build_circuitry_polytope(
    tensor = tensor,
    embedding = embedding,
    circuitry_id = circuitry_id
  )
  
  plot_circuitry_polytope_custom(
    poly_obj = poly,
    title = paste0("Panel ", panel_label, " — ", circuitry_id),
    outfile_html = NULL
  )
}

fig_list_4 <- setNames(vector("list", 4), panel_order)
for (k in seq_len(nrow(tbl))) {
  lbl <- tbl$panel[k]
  cid <- tbl$Circuitries_id[k]
  fig_list_4[[lbl]] <- build_panel_plotly(lbl, cid, tensor, embedding)
}

# ---------------------------------------------------------------------------
# 2) Latent dimension names (18D) — shown ONCE as a top-right block
# ---------------------------------------------------------------------------

latent_dim_names <- c(
  "rho", "rho_strength", "TN_dir", "TN_strength",
  "OS_dir", "OS_strength", "OS_lr_chisq",
  "DSS_dir", "DSS_strength", "DSS_lr_chisq",
  "DFI_dir", "DFI_strength", "DFI_lr_chisq",
  "PFI_dir", "PFI_strength", "PFI_lr_chisq",
  "TME_score", "Immune_dir"
)
stopifnot(length(latent_dim_names) == 18)

# ---------------------------------------------------------------------------
# 3) TRUE 2×2 multiscene composer
#    - Each panel goes to its own scene: scene, scene2, scene3, scene4
#    - Legend appears ONCE (from first panel only)
#    - Right margin reserved for the two-tier legend area
# ---------------------------------------------------------------------------

make_multiscene_3d <- function(fig_list,
                               nrows = 2, ncols = 2,
                               margin_x = 0.04, margin_y = 0.08,
                               latent_dim_names = NULL,
                               right_margin_px = 360,
                               top_margin_px = 70,
                               bottom_margin_px = 10,
                               left_margin_px = 10,
                               # Right-side placement (paper coordinates)
                               legend_x = 1.02,
                               legend_y = 0.40,      # trace legend: LOWER block
                               dimblock_x = 1.02,
                               dimblock_y = 0.98,    # 18D block: TOP block (must be <= 1)
                               dimblock_fontsize = 11,
                               legend_fontsize = 11) {
  
  stopifnot(length(fig_list) == nrows * ncols)
  
  # ---- domains for each scene ----------------------------------------------
  x_breaks <- seq(0, 1, length.out = ncols + 1)
  y_breaks <- seq(0, 1, length.out = nrows + 1)
  
  x_domains <- lapply(seq_len(ncols), function(j) {
    c(x_breaks[j] + margin_x/2, x_breaks[j+1] - margin_x/2)
  })
  
  y_domains <- lapply(seq_len(nrows), function(i) {
    top_row_index <- nrows - i + 1
    c(y_breaks[top_row_index] + margin_y/2, y_breaks[top_row_index+1] - margin_y/2)
  })
  
  # ---- inject traces into distinct scenes -----------------------------------
  # <<< PATCH (optional but safe): initialize as a 3D container to avoid
  # Plotly creating a default 2D cartesian subplot with giant X/Y axes.
  out <- plotly::plot_ly(type = "scatter3d")  # <<< PATCH
  out$x$data <- list()
  
  for (k in seq_along(fig_list)) {
    
    pb <- plotly::plotly_build(fig_list[[k]])
    scene_id <- if (k == 1) "scene" else paste0("scene", k)
    
    for (tr in pb$x$data) {
      
      # assign scene
      tr$scene <- scene_id
      
      # IMPORTANT: show the plotly legend only once (from panel 1)
      # This prevents a 4× duplicated legend and also stops overlap pressure.
      if (k != 1) {
        tr$showlegend <- FALSE
      }
      
      out$x$data <- c(out$x$data, list(tr))
    }
  }
  
  # ---- base layout ----------------------------------------------------------
  lay <- list(
    margin = list(
      l = left_margin_px,
      r = right_margin_px,
      b = bottom_margin_px,
      t = top_margin_px
    ),
    
    # <<< PATCH (critical): disable any global 2D cartesian axes that Plotly
    # may attach to the root figure (the “giant” X/Y you saw).
    xaxis = list(visible = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(visible = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    
    # Single trace legend (LOWER right block)
    legend = list(
      x = legend_x,
      y = legend_y,
      xref = "paper",
      yref = "paper",
      xanchor = "left",
      yanchor = "top",
      orientation = "v",
      bgcolor = "rgba(255,255,255,0.88)",
      bordercolor = "rgba(0,0,0,0.15)",
      borderwidth = 1,
      font = list(size = legend_fontsize)
    )
  )
  
  # ---- per-scene layout + per-panel titles ---------------------------------
  for (k in seq_along(fig_list)) {
    
    scene_id <- if (k == 1) "scene" else paste0("scene", k)
    
    row_i <- ceiling(k / ncols)
    col_j <- k - (row_i - 1) * ncols
    
    pb <- plotly::plotly_build(fig_list[[k]])
    
    this_title <- ""
    if (!is.null(pb$x$layout$title) && !is.null(pb$x$layout$title$text)) {
      this_title <- pb$x$layout$title$text
    }
    
    lay[[scene_id]] <- list(
      domain = list(x = x_domains[[col_j]], y = y_domains[[row_i]]),
      xaxis  = list(title = "PC1"),
      yaxis  = list(title = "PC2"),
      zaxis  = list(title = "PC3")
    )
    
    # panel title above each scene (annotation)
    if (nzchar(this_title)) {
      if (is.null(lay$annotations)) lay$annotations <- list()
      lay$annotations <- c(lay$annotations, list(
        list(
          text = this_title,
          x = mean(x_domains[[col_j]]),
          y = y_domains[[row_i]][2] + 0.015,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          font = list(size = 12),
          xanchor = "center",
          yanchor = "bottom",
          align = "center"
        )
      ))
    }
  }
  
  # ---- 18D latent axes block (TOP right block) ------------------------------
  if (!is.null(latent_dim_names)) {
    
    stopifnot(length(latent_dim_names) == 18)
    
    dim_txt <- paste(
      "<b>Latent dimensions (18D):</b>",
      paste(sprintf("%02d: %s", seq_along(latent_dim_names), latent_dim_names),
            collapse = "<br>"),
      sep = "<br>"
    )
    
    if (is.null(lay$annotations)) lay$annotations <- list()
    lay$annotations <- c(lay$annotations, list(
      list(
        text = dim_txt,
        x = dimblock_x,
        y = dimblock_y,
        xref = "paper",
        yref = "paper",
        xanchor = "left",
        yanchor = "top",
        align = "left",
        showarrow = FALSE,
        bgcolor = "rgba(255,255,255,0.92)",
        bordercolor = "rgba(0,0,0,0.15)",
        borderwidth = 1,
        font = list(size = dimblock_fontsize)
      )
    ))
  }
  
  # ---- apply layout safely --------------------------------------------------
  out <- do.call(plotly::layout, c(list(out), lay))
  out
}

# ---------------------------------------------------------------------------
# 4) Build composite + save + open (robust on Windows)
# ---------------------------------------------------------------------------

fig4_multi <- make_multiscene_3d(
  fig_list = list(fig_list_4$A, fig_list_4$B, fig_list_4$C, fig_list_4$D),
  nrows = 2, ncols = 2,
  latent_dim_names = latent_dim_names,
  # These two numbers are the “no-overlap” controls:
  dimblock_y = 0.98,   # TOP block (<= 1)
  legend_y   = 0.40,   # LOWER block
  right_margin_px = 380,
  top_margin_px = 70
)

out_dir <- "Figure_4_panels_html"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "Figure_4_MULTISCENE_2x2.html")
htmlwidgets::saveWidget(fig4_multi, out_file, selfcontained = TRUE)

# Robust open on Windows (important: file:/// + forward slashes)
browseURL(paste0("file:///", normalizePath(out_file, winslash = "/")))

####
####
#### Figure S2
####
####
####

library(dplyr)
library(stringr)

# -------------------------------------------------------------------------
# Helper predicates for your *actual* classifier encoding
# -------------------------------------------------------------------------
is_discordant <- function(x) {
  x %in% c("extreme_discordance", "strong_discordance", "moderate_discordance")
}

is_low_dim <- function(x) {
  str_detect(x, "^low_dimensional_flat\\s*;")
}

is_intermediate <- function(x) {
  str_detect(x, "^intermediate_complexity\\s*;")
}

is_high_dim <- function(x) {
  str_detect(x, "^high_complexity_multidimensional\\s*;")
}

# Optional: restrict to strongly asymmetric only (if you want "clean" exemplars)
is_strong_asym <- function(x) {
  str_detect(x, "strongly_asymmetric_sig_int")
}

# -------------------------------------------------------------------------
# Deterministic representative selector: closest to the median barycenter_distance
# within a candidate set (prevents cherry-picking and avoids random sampling)
# -------------------------------------------------------------------------
select_median_rep <- function(df) {
  stopifnot("barycenter_distance" %in% names(df))
  med <- median(df$barycenter_distance, na.rm = TRUE)
  
  df %>%
    mutate(.delta = abs(barycenter_distance - med)) %>%
    arrange(.delta) %>%
    slice(1) %>%
    select(-.delta)
}

# -------------------------------------------------------------------------
# Figure S2 panels (A–D): ONE circuitry per panel
# NOTE: This implements the population-level taxonomy logic.
# -------------------------------------------------------------------------

# S2A: Discordant + low-dimensional (geometrically simple, separated barycenters)
S2A <- circuitries %>%
  filter(is_discordant(distance_implication), is_low_dim(vol_implication)) %>%
  select_median_rep()

# S2B: Concordant + low-dimensional (compact, overlapping)
S2B <- circuitries %>%
  filter(distance_implication == "high_concordance", is_low_dim(vol_implication)) %>%
  select_median_rep()

# S2C: Discordant + intermediate-complexity (often component-dominated)
S2C <- circuitries %>%
  filter(is_discordant(distance_implication), is_intermediate(vol_implication)) %>%
  select_median_rep()

# S2D: Discordant + high-complexity multidimensional (dominant atlas phenotype)
S2D <- circuitries %>%
  filter(is_discordant(distance_implication), is_high_dim(vol_implication)) %>%
  select_median_rep()

# -------------------------------------------------------------------------
# Combine into a single data frame for export / inspection
# -------------------------------------------------------------------------
S2_reps <- bind_rows(
  S2A %>% mutate(FigureS2_panel = "A"),
  S2B %>% mutate(FigureS2_panel = "B"),
  S2C %>% mutate(FigureS2_panel = "C"),
  S2D %>% mutate(FigureS2_panel = "D")
)

# Quick sanity check output: show the selected IDs and regimes
S2_reps %>%
  select(FigureS2_panel, Circuitries_id, distance_implication, vol_implication, barycenter_distance)


####
####
#### Checking the selected Figure S2 panels individually
####
####

####
# # Or by specific Circuitries_id:
poly_ESCA <- build_circuitry_polytope(
  tensor,
  embedding,
  circuitry_id = "ESCA-2471 / ESCA-4395"
)
plot_circuitry_polytope(poly_ESCA)

# # Or by specific Circuitries_id:
poly_KICH <- build_circuitry_polytope(
  tensor,
  embedding,
  circuitry_id = "KICH-3614 / KICH-3961"
)
plot_circuitry_polytope(poly_KICH)

# # Or by specific Circuitries_id:
poly_BLCA <- build_circuitry_polytope(
  tensor,
  embedding,
  circuitry_id = "BLCA-4374 / BLCA-7378"
)
plot_circuitry_polytope(poly_BLCA)

# # Or by specific Circuitries_id:
poly_PAAD <- build_circuitry_polytope(
  tensor,
  embedding,
  circuitry_id = "PAAD-5314 / PAAD-2498"
)
plot_circuitry_polytope(poly_PAAD)





###############################################################################
# Figure S2 Final — TRUE 2×2 multiscene HTML
# FIX: Two-line per-panel headings to avoid vertical collisions
#   Line 1 (bold):  "Panel A — ESCA-2471 / ESCA-4395"
#   Line 2 (smaller): "strong_discordance | low_dimensional_flat ; ..."
###############################################################################

stopifnot(requireNamespace("plotly", quietly = TRUE))
stopifnot(requireNamespace("htmlwidgets", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("stringr", quietly = TRUE))

library(dplyr)
library(stringr)

# ---- Required objects -------------------------------------------------------
stopifnot(exists("S2_reps"))
stopifnot(all(c("Circuitries_id","distance_implication","vol_implication") %in% names(S2_reps)))
stopifnot(exists("tensor"))
stopifnot(exists("embedding"))
stopifnot(exists("build_circuitry_polytope"))
stopifnot(exists("plot_circuitry_polytope_custom"))

# ---------------------------------------------------------------------------
# 0) Output directory
# ---------------------------------------------------------------------------
out_dir_S2 <- "Figure_S2_panels_html"
dir.create(out_dir_S2, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1) Tag S2 rows to panels A–D by regime (NOT by row order)
# ---------------------------------------------------------------------------
is_low_dim <- function(x) stringr::str_detect(x, "^low_dimensional_flat\\s*;")
is_inter   <- function(x) stringr::str_detect(x, "^intermediate_complexity\\s*;")
is_high    <- function(x) stringr::str_detect(x, "^high_complexity_multidimensional\\s*;")
is_discord <- function(x) x %in% c("extreme_discordance","strong_discordance","moderate_discordance")

S2_tagged <- S2_reps %>%
  mutate(
    panel = dplyr::case_when(
      is_discord(distance_implication) & is_low_dim(vol_implication) ~ "A",
      distance_implication == "high_concordance" & is_low_dim(vol_implication) ~ "B",
      is_discord(distance_implication) & is_inter(vol_implication) ~ "C",
      is_discord(distance_implication) & is_high(vol_implication) ~ "D",
      TRUE ~ NA_character_
    )
  )

stopifnot(all(!is.na(S2_tagged$panel)))
stopifnot(n_distinct(S2_tagged$panel) == 4)

make_key <- function(x) {
  x %>%
    stringr::str_replace_all("\\s*/\\s*", "_") %>%
    stringr::str_replace_all("\\s+", "_") %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]", "")
}

report_S2 <- S2_tagged %>%
  transmute(
    panel = as.character(panel),
    Circuitries_id = as.character(Circuitries_id),
    circuitry_key = make_key(Circuitries_id),
    html_file = file.path(out_dir_S2, paste0("Figure_S2_Panel_", panel, "_", circuitry_key, ".html")),
    distance_implication = as.character(distance_implication),
    vol_implication = as.character(vol_implication)
  )

panel_order <- c("A","B","C","D")
report_S2 <- report_S2[order(match(report_S2$panel, panel_order)), ]
stopifnot(nrow(report_S2) == 4)

# ---------------------------------------------------------------------------
# 2) Two-line panel titles (collision-safe)
#    - Line 1: Panel X — Circuitries_id
#    - Line 2: distance_implication | vol_implication   (smaller font)
# ---------------------------------------------------------------------------
make_S2_panel_title <- function(panel_label, circuitry_id, dist_imp, vol_imp) {
  line1 <- sprintf("Panel %s — %s", panel_label, circuitry_id)
  
  # Split line 2 into TWO short lines: distance on one line, volume on the next.
  # This avoids any midline collision and does not depend on viewport width.
  line2a <- sprintf("%s", dist_imp)
  line2b <- sprintf("%s", vol_imp)
  
  paste0(
    "<b>", line1, "</b>",
    "<br><span style='font-size:10px'>", line2a, "</span>",
    "<br><span style='font-size:10px'>", line2b, "</span>"
  )
}

panel_titles_S2 <- vector("list", 4)
for (i in seq_len(nrow(report_S2))) {
  panel_titles_S2[[i]] <- make_S2_panel_title(
    panel_label = report_S2$panel[i],
    circuitry_id = report_S2$Circuitries_id[i],
    dist_imp = report_S2$distance_implication[i],
    vol_imp  = report_S2$vol_implication[i]
  )
}

# ---------------------------------------------------------------------------
# 3) Latent dimension names (18D) — same as Figure 4
# ---------------------------------------------------------------------------
latent_dim_names <- c(
  "rho", "rho_strength", "TN_dir", "TN_strength",
  "OS_dir", "OS_strength", "OS_lr_chisq",
  "DSS_dir", "DSS_strength", "DSS_lr_chisq",
  "DFI_dir", "DFI_strength", "DFI_lr_chisq",
  "PFI_dir", "PFI_strength", "PFI_lr_chisq",
  "TME_score", "Immune_dir"
)
stopifnot(length(latent_dim_names) == 18)

# ---------------------------------------------------------------------------
# 4) Build the FOUR plotly objects in memory (A–D)
#    NOTE: We keep the plot's internal title empty; we supply our own annotations.
# ---------------------------------------------------------------------------
build_panel_plotly_S2 <- function(circuitry_id, tensor, embedding) {
  
  poly <- build_circuitry_polytope(
    tensor = tensor,
    embedding = embedding,
    circuitry_id = circuitry_id
  )
  
  # IMPORTANT: title = "" so we don't get competing titles from pb$x$layout$title$text
  plot_circuitry_polytope_custom(
    poly_obj = poly,
    title = "",
    outfile_html = NULL
  )
}

fig_list_S2 <- setNames(vector("list", 4), panel_order)
for (k in seq_len(nrow(report_S2))) {
  lbl <- report_S2$panel[k]
  cid <- report_S2$Circuitries_id[k]
  fig_list_S2[[lbl]] <- build_panel_plotly_S2(cid, tensor, embedding)
}

# ---------------------------------------------------------------------------
# 5) TRUE 2×2 multiscene composer (Figure 4 engine + title injection)
#    - Legends as before (18D block + one trace legend)
#    - NEW: panel_titles (two-line strings) and smaller y-offset
# ---------------------------------------------------------------------------
make_multiscene_3d <- function(fig_list,
                               nrows = 2, ncols = 2,
                               margin_x = 0.08, margin_y = 0.08,
                               latent_dim_names = NULL,
                               right_margin_px = 470,
                               top_margin_px = 70,
                               bottom_margin_px = 10,
                               left_margin_px = 10,
                               legend_x = 1.02,
                               legend_y = 0.40,
                               dimblock_x = 1.02,
                               dimblock_y = 0.98,
                               dimblock_fontsize = 11,
                               legend_fontsize = 11,
                               panel_titles = NULL,
                               panel_title_y_offset = 0.006) {
  
  stopifnot(length(fig_list) == nrows * ncols)
  if (!is.null(panel_titles)) stopifnot(length(panel_titles) == length(fig_list))
  
  # ---- domains for each scene ----------------------------------------------
  x_breaks <- seq(0, 1, length.out = ncols + 1)
  y_breaks <- seq(0, 1, length.out = nrows + 1)
  
  x_domains <- lapply(seq_len(ncols), function(j) {
    c(x_breaks[j] + margin_x/2, x_breaks[j+1] - margin_x/2)
  })
  
  y_domains <- lapply(seq_len(nrows), function(i) {
    top_row_index <- nrows - i + 1
    c(y_breaks[top_row_index] + margin_y/2, y_breaks[top_row_index+1] - margin_y/2)
  })
  
  # ---- inject traces into distinct scenes -----------------------------------
  out <- plotly::plot_ly(type = "scatter3d")  # safe 3D container
  out$x$data <- list()
  
  for (k in seq_along(fig_list)) {
    
    pb <- plotly::plotly_build(fig_list[[k]])
    scene_id <- if (k == 1) "scene" else paste0("scene", k)
    
    for (tr in pb$x$data) {
      
      tr$scene <- scene_id
      
      # show trace legend only once (panel 1)
      if (k != 1) tr$showlegend <- FALSE
      
      out$x$data <- c(out$x$data, list(tr))
    }
  }
  
  # ---- base layout ----------------------------------------------------------
  lay <- list(
    margin = list(
      l = left_margin_px,
      r = right_margin_px,
      b = bottom_margin_px,
      t = top_margin_px
    ),
    
    # Disable any global 2D cartesian axes Plotly might attach
    xaxis = list(visible = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(visible = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    
    legend = list(
      x = legend_x, y = legend_y, xref = "paper", yref = "paper",
      xanchor = "left", yanchor = "top",
      orientation = "v",
      bgcolor = "rgba(255,255,255,0.88)",
      bordercolor = "rgba(0,0,0,0.15)",
      borderwidth = 1,
      font = list(size = legend_fontsize)
    )
  )
  
  # ---- per-scene layout + custom per-panel titles ---------------------------
  for (k in seq_along(fig_list)) {
    
    scene_id <- if (k == 1) "scene" else paste0("scene", k)
    
    row_i <- ceiling(k / ncols)
    col_j <- k - (row_i - 1) * ncols
    
    lay[[scene_id]] <- list(
      domain = list(x = x_domains[[col_j]], y = y_domains[[row_i]]),
      xaxis  = list(title = "PC1"),
      yaxis  = list(title = "PC2"),
      zaxis  = list(title = "PC3")
    )
    
    # Use provided two-line title (HTML)
    if (!is.null(panel_titles)) {
      if (is.null(lay$annotations)) lay$annotations <- list()
      lay$annotations <- c(lay$annotations, list(
        list(
          text = panel_titles[[k]],
          x = mean(x_domains[[col_j]]),
          y = y_domains[[row_i]][2] + panel_title_y_offset,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          font = list(size = 12),
          xanchor = "center",
          yanchor = "bottom",
          align = "center"
        )
      ))
    }
  }
  
  # ---- 18D latent axes block (TOP right block) ------------------------------
  if (!is.null(latent_dim_names)) {
    
    stopifnot(length(latent_dim_names) == 18)
    
    dim_txt <- paste(
      "<b>Latent dimensions (18D):</b>",
      paste(sprintf("%02d: %s", seq_along(latent_dim_names), latent_dim_names),
            collapse = "<br>"),
      sep = "<br>"
    )
    
    if (is.null(lay$annotations)) lay$annotations <- list()
    lay$annotations <- c(lay$annotations, list(
      list(
        text = dim_txt,
        x = dimblock_x,
        y = dimblock_y,
        xref = "paper",
        yref = "paper",
        xanchor = "left",
        yanchor = "top",
        align = "left",
        showarrow = FALSE,
        bgcolor = "rgba(255,255,255,0.92)",
        bordercolor = "rgba(0,0,0,0.15)",
        borderwidth = 1,
        font = list(size = dimblock_fontsize)
      )
    ))
  }
  
  out <- do.call(plotly::layout, c(list(out), lay))
  out
}

# ---------------------------------------------------------------------------
# 6) Build composite + save + open
# ---------------------------------------------------------------------------
figS2_multi <- make_multiscene_3d(
  fig_list = list(fig_list_S2$A, fig_list_S2$B, fig_list_S2$C, fig_list_S2$D),
  nrows = 2, ncols = 2,
  
  margin_x = 0.08,   # << was 0.04; increases horizontal gap between columns
  margin_y = 0.08,   # keep as-is unless you also want more row gap
  
  latent_dim_names = latent_dim_names,
  
  # move legend column further right (see Patch 2)
  right_margin_px = 470,
  legend_x = 1.08,
  dimblock_x = 1.08,
  
  dimblock_y = 0.98,
  legend_y   = 0.40,
  
  top_margin_px = 70,
  panel_titles = panel_titles_S2,
  panel_title_y_offset = 0.006
)

out_file_S2 <- file.path(out_dir_S2, "Figure_S2_MULTISCENE_2x2.html")
htmlwidgets::saveWidget(figS2_multi, out_file_S2, selfcontained = TRUE)

browseURL(paste0("file:///", normalizePath(out_file_S2, winslash = "/")))

webshot2::webshot(
  "Figure_S2_panels_html/Figure_S2_MULTISCENE_2x2.html",
  file = "Figure_S2_MULTISCENE_2x2.png",
  vwidth = 2400,
  vheight = 3000,
  zoom = 2.5
)

magick::image_read("Figure_S2_MULTISCENE_2x2.png") |>
  magick::image_write("Figure_S2_MULTISCENE_2x2_600dpi.tiff", density = "600x600")

# Convert PNG → 600 dpi TIFF
magick::image_write(
  magick::image_read("Figure_S2_MULTISCENE_2x2.png"),
  path    = "Figure_S2_MULTISCENE_2x2_600dpi.tiff",
  format  = "tiff",
  density = "600x600"
)


#####
#####
#### Extracting relevant variables and values for Figure S2 multipanel legend
####
####
####
###############################################################################
# EXTRACT RELEVANT VARIABLES FOR THE CIRCUITRY: ESCA-2471 / ESCA-4395
###############################################################################

# 1. Target circuitry identifier
target_id <- "ESCA-2471 / ESCA-4395"

# 2. Variables to extract (same structure as cir_relevant_variables.tsv)
vars_to_extract <- c(
  "Circuitries_id",
  "Nomenclature_sig",
  "Signatures",
  "Nomenclature_int",
  "Interaction",
  "barycenter_distance",
  "sig_hull_vol",
  "int_hull_vol",
  "vol_ratio",
  "distance_implication",
  "vol_implication",
  "Omic_layer_sig",
  "Phenotypic_layer_sig",
  "Omic_layer_int",
  "Phenotypic_layer_int",
  "Correlation_rho_sig",
  "Correlation_rho_int",
  "Tumor_vs_normal_sig",
  "Tumor_vs_normal_int",
  "Combined_outcome_HRC_sig",
  "Combined_outcome_HRC_int",
  "Microenvironment_classification_sig",
  "Microenvironment_classification_int",
  "Immune_classification_sig",
  "Immune_classification_int",
  "Phenotypic_concordance",
  "Survival_concordance_aggregated",
  "Immune_concordance",
  "Final_concordance_summary",
  "CTAB",
  "Metabolism",
  "Pathways"
)
# 3. Extract, convert all values to character (required for tidy pivot), pivot to two-column format
cir_2471_4395_relevant_variables <- circuitries %>%
  dplyr::filter(Circuitries_id == target_id) %>%
  dplyr::select(dplyr::all_of(vars_to_extract)) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%  # ensure homogeneous type
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# 4. Print result
cir_2471_4395_relevant_variables

# 5. (Optional) Save to TSV
readr::write_tsv(cir_2471_4395_relevant_variables, "ESCA_2471_ESCA_4395_cir_relevant_variables.tsv")
###############################################################################

###############################################################################
# EXTRACT RELEVANT VARIABLES FOR THE CIRCUITRY: KICH-3614 / KICH-3961
###############################################################################

# 1. Target circuitry identifier
target_id <- "KICH-3614 / KICH-3961"

# 2. Variables to extract (same structure as cir_relevant_variables.tsv)
vars_to_extract <- c(
  "Circuitries_id",
  "Nomenclature_sig",
  "Signatures",
  "Nomenclature_int",
  "Interaction",
  "barycenter_distance",
  "sig_hull_vol",
  "int_hull_vol",
  "vol_ratio",
  "distance_implication",
  "vol_implication",
  "Omic_layer_sig",
  "Phenotypic_layer_sig",
  "Omic_layer_int",
  "Phenotypic_layer_int",
  "Correlation_rho_sig",
  "Correlation_rho_int",
  "Tumor_vs_normal_sig",
  "Tumor_vs_normal_int",
  "Combined_outcome_HRC_sig",
  "Combined_outcome_HRC_int",
  "Microenvironment_classification_sig",
  "Microenvironment_classification_int",
  "Immune_classification_sig",
  "Immune_classification_int",
  "Phenotypic_concordance",
  "Survival_concordance_aggregated",
  "Immune_concordance",
  "Final_concordance_summary",
  "CTAB",
  "Metabolism",
  "Pathways"
)
# 3. Extract, convert all values to character (required for tidy pivot), pivot to two-column format
cir_3614_3961_relevant_variables <- circuitries %>%
  dplyr::filter(Circuitries_id == target_id) %>%
  dplyr::select(dplyr::all_of(vars_to_extract)) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%  # ensure homogeneous type
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# 4. Print result
cir_3614_3961_relevant_variables

# 5. (Optional) Save to TSV
readr::write_tsv(cir_3614_3961_relevant_variables, "KICH_3614_KICH_3961_cir_relevant_variables.tsv")
###############################################################################

###############################################################################
# EXTRACT RELEVANT VARIABLES FOR THE CIRCUITRY: BLCA-4374 / BLCA-7378
###############################################################################

# 1. Target circuitry identifier
target_id <- "BLCA-4374 / BLCA-7378"

# 2. Variables to extract (same structure as cir_relevant_variables.tsv)
vars_to_extract <- c(
  "Circuitries_id",
  "Nomenclature_sig",
  "Signatures",
  "Nomenclature_int",
  "Interaction",
  "barycenter_distance",
  "sig_hull_vol",
  "int_hull_vol",
  "vol_ratio",
  "distance_implication",
  "vol_implication",
  "Omic_layer_sig",
  "Phenotypic_layer_sig",
  "Omic_layer_int",
  "Phenotypic_layer_int",
  "Correlation_rho_sig",
  "Correlation_rho_int",
  "Tumor_vs_normal_sig",
  "Tumor_vs_normal_int",
  "Combined_outcome_HRC_sig",
  "Combined_outcome_HRC_int",
  "Microenvironment_classification_sig",
  "Microenvironment_classification_int",
  "Immune_classification_sig",
  "Immune_classification_int",
  "Phenotypic_concordance",
  "Survival_concordance_aggregated",
  "Immune_concordance",
  "Final_concordance_summary",
  "CTAB",
  "Metabolism",
  "Pathways"
)
# 3. Extract, convert all values to character (required for tidy pivot), pivot to two-column format
cir_4374_73781_relevant_variables <- circuitries %>%
  dplyr::filter(Circuitries_id == target_id) %>%
  dplyr::select(dplyr::all_of(vars_to_extract)) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%  # ensure homogeneous type
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# 4. Print result
cir_4374_73781_relevant_variables

# 5. (Optional) Save to TSV
readr::write_tsv(cir_4374_73781_relevant_variables, "BLCA_4374_BKCA_7378_cir_relevant_variables.tsv")
###############################################################################
###############################################################################
# EXTRACT RELEVANT VARIABLES FOR THE CIRCUITRY: PAAD-5314 / PAAD-2498
###############################################################################

# 1. Target circuitry identifier
target_id <- "PAAD-5314 / PAAD-2498"

# 2. Variables to extract (same structure as cir_relevant_variables.tsv)
vars_to_extract <- c(
  "Circuitries_id",
  "Nomenclature_sig",
  "Signatures",
  "Nomenclature_int",
  "Interaction",
  "barycenter_distance",
  "sig_hull_vol",
  "int_hull_vol",
  "vol_ratio",
  "distance_implication",
  "vol_implication",
  "Omic_layer_sig",
  "Phenotypic_layer_sig",
  "Omic_layer_int",
  "Phenotypic_layer_int",
  "Correlation_rho_sig",
  "Correlation_rho_int",
  "Tumor_vs_normal_sig",
  "Tumor_vs_normal_int",
  "Combined_outcome_HRC_sig",
  "Combined_outcome_HRC_int",
  "Microenvironment_classification_sig",
  "Microenvironment_classification_int",
  "Immune_classification_sig",
  "Immune_classification_int",
  "Phenotypic_concordance",
  "Survival_concordance_aggregated",
  "Immune_concordance",
  "Final_concordance_summary",
  "CTAB",
  "Metabolism",
  "Pathways"
)
# 3. Extract, convert all values to character (required for tidy pivot), pivot to two-column format
cir_5314_2498_relevant_variables <- circuitries %>%
  dplyr::filter(Circuitries_id == target_id) %>%
  dplyr::select(dplyr::all_of(vars_to_extract)) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%  # ensure homogeneous type
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Value"
  )

# 4. Print result
cir_5314_2498_relevant_variables

# 5. (Optional) Save to TSV
readr::write_tsv(cir_5314_2498_relevant_variables, "PAAD_5314_PAAD_2498_cir_relevant_variables.tsv")
###############################################################################

#####
#####
##### Figure S1
#####
#####
#####

###############################################################################
# Figure S1 — Select 6 representative circuitries (3×2 classes) using
# vol_implication directly (complexity × symmetry-state)
#
# FIX IMPLEMENTED:
# - Avoid substring collision where "symmetric_sig_int" is contained within
#   "strongly_asymmetric_sig_int" and "moderately_asymmetric_sig_int".
# - Parse vol_implication by splitting on ';' and matching tokens EXACTLY.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

select_fs1_representatives <- function(circuitries) {
  
  # ---- Input checks ----
  req <- c("Circuitries_id", "sig_hull_vol", "int_hull_vol", "vol_ratio", "vol_implication")
  missing_cols <- setdiff(req, names(circuitries))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # ---- Normalize and split vol_implication into exact tokens ----
  parsed <- circuitries %>%
    mutate(
      vol_implication_clean = str_squish(vol_implication),
      vol_tokens = str_split(vol_implication_clean, "\\s*;\\s*")
    ) %>%
    mutate(
      # Extract tokens robustly (exact values)
      token_complexity = vapply(vol_tokens, function(x) if (length(x) >= 1) x[[1]] else NA_character_, character(1)),
      token_asymmetry  = vapply(vol_tokens, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
    ) %>%
    mutate(
      # Complexity tier from token 1 (exact match)
      complexity_tier = case_when(
        token_complexity == "low_dimensional_flat" ~ "Low",
        token_complexity == "intermediate_complexity" ~ "Intermediate",
        token_complexity == "high_complexity_multidimensional" ~ "High",
        TRUE ~ NA_character_
      ),
      
      # Symmetry-state (2-level for Figure S1) from token 2 (exact match)
      symmetry_state = case_when(
        token_asymmetry == "symmetric_sig_int" ~ "Symmetric",
        token_asymmetry %in% c("moderately_asymmetric_sig_int", "strongly_asymmetric_sig_int") ~ "Sig_dominant",
        TRUE ~ NA_character_
      ),
      
      class_3x2 = if_else(
        !is.na(complexity_tier) & !is.na(symmetry_state),
        paste0(complexity_tier, "_", symmetry_state),
        NA_character_
      )
    )
  
  # ---- Validate parsing completeness ----
  n_na <- sum(is.na(parsed$class_3x2))
  if (n_na > 0) {
    warning(
      "Some rows could not be assigned to a 3×2 class from vol_implication (n = ",
      n_na, "). They will be excluded from selection."
    )
  }
  
  parsed2 <- parsed %>%
    filter(!is.na(class_3x2)) %>%
    filter(
      is.finite(sig_hull_vol), is.finite(int_hull_vol), is.finite(vol_ratio),
      sig_hull_vol > 0, int_hull_vol > 0, vol_ratio > 0
    )
  
  # ---- Log-space centrality ----
  eps <- 1e-12
  parsed2 <- parsed2 %>%
    mutate(
      log_sig_vol   = log10(sig_hull_vol + eps),
      log_int_vol   = log10(int_hull_vol + eps),
      log_vol_ratio = log10(vol_ratio + eps)
    )
  
  # ---- Select class-central exemplar per 3×2 class ----
  representatives_df <- parsed2 %>%
    group_by(class_3x2) %>%
    mutate(
      med_log_sig    = median(log_sig_vol, na.rm = TRUE),
      med_log_int    = median(log_int_vol, na.rm = TRUE),
      med_log_ratio  = median(log_vol_ratio, na.rm = TRUE),
      dist_to_median = sqrt(
        (log_sig_vol   - med_log_sig)^2 +
          (log_int_vol   - med_log_int)^2 +
          (log_vol_ratio - med_log_ratio)^2
      )
    ) %>%
    arrange(dist_to_median) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      grid_row = factor(complexity_tier, levels = c("Low", "Intermediate", "High")),
      grid_col = factor(symmetry_state, levels = c("Symmetric", "Sig_dominant"))
    ) %>%
    select(
      class_3x2, complexity_tier, symmetry_state, grid_row, grid_col,
      Circuitries_id, sig_hull_vol, int_hull_vol, vol_ratio, vol_implication,
      dist_to_median
    ) %>%
    arrange(grid_row, grid_col)
  
  # ---- Check expected 6 classes ----
  expected_classes <- as.vector(outer(
    c("Low", "Intermediate", "High"),
    c("Symmetric", "Sig_dominant"),
    paste, sep = "_"
  ))
  missing_classes <- setdiff(expected_classes, representatives_df$class_3x2)
  if (length(missing_classes) > 0) {
    warning(
      "Not all 6 (3×2) classes were found in the data. Missing classes: ",
      paste(missing_classes, collapse = ", ")
    )
  }
  
  return(representatives_df)
}

###############################################################################
# USAGE:
representatives_df <- select_fs1_representatives(circuitries)
print(dim(representatives_df))   # should now be 6 × 11
print(representatives_df)
###############################################################################

representatives_df <- representatives_df %>%
  mutate(
    panel = c("A","B","C","D","E","F"),
    panel_label = paste0(
      panel, ": ",
      as.character(grid_row), " — ",
      ifelse(as.character(grid_col) == "Symmetric", "Symmetric", "Signature-dominant")
    )
  )

representatives_df %>%
  select(panel, panel_label, class_3x2, Circuitries_id, vol_implication,
         sig_hull_vol, int_hull_vol, vol_ratio) %>%
  print(n = Inf)

rio::export(representatives_df, "Figure_S1_representative_df.tsv")

###############################################################################
# Supplementary Table S1 — Figure S1 Panel Key
# Builds a clean panel-mapping table from representatives_df (6 × 13)
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

build_fs1_panel_key <- function(representatives_df) {
  
  req <- c("panel","panel_label","class_3x2","complexity_tier","symmetry_state",
           "Circuitries_id","sig_hull_vol","int_hull_vol","vol_ratio",
           "vol_implication","dist_to_median")
  missing_cols <- setdiff(req, names(representatives_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in representatives_df: ", paste(missing_cols, collapse = ", "))
  }
  
  out <- representatives_df %>%
    mutate(
      # Parse Circuitries_id: e.g., "PAAD-5314 / PAAD-2498"
      Circuitry_sig_id = str_trim(str_split_fixed(Circuitries_id, "/", 2)[,1]),
      Circuitry_int_id = str_trim(str_split_fixed(Circuitries_id, "/", 2)[,2]),
      CTAB = str_extract(Circuitry_sig_id, "^[A-Za-z0-9]+"),
      # Pretty class label for human readers
      Geometry_class = case_when(
        class_3x2 == "Low_Symmetric" ~ "Low-dimensional, symmetric geometry",
        class_3x2 == "Low_Sig_dominant" ~ "Low-dimensional, signature-dominant asymmetric geometry",
        class_3x2 == "Intermediate_Symmetric" ~ "Intermediate-complexity, symmetric geometry",
        class_3x2 == "Intermediate_Sig_dominant" ~ "Intermediate-complexity, signature-dominant asymmetric geometry",
        class_3x2 == "High_Symmetric" ~ "High-complexity multidimensional, symmetric geometry",
        class_3x2 == "High_Sig_dominant" ~ "High-complexity multidimensional, signature-dominant asymmetric geometry",
        TRUE ~ class_3x2
      )
    ) %>%
    # Format numeric columns (keep raw + formatted if you want; here we format for table)
    mutate(
      sig_hull_vol_fmt = format(sig_hull_vol, scientific = TRUE, digits = 4),
      int_hull_vol_fmt = format(int_hull_vol, scientific = TRUE, digits = 4),
      vol_ratio_fmt    = format(vol_ratio, scientific = TRUE, digits = 4),
      dist_to_median_fmt = format(dist_to_median, scientific = TRUE, digits = 4)
    ) %>%
    transmute(
      Panel = panel,
      Panel_description = panel_label,
      CTAB = CTAB,
      Circuitries_id = Circuitries_id,
      Circuitry_sig_id = Circuitry_sig_id,
      Circuitry_int_id = Circuitry_int_id,
      Geometry_class = Geometry_class,
      Complexity_tier = complexity_tier,
      Symmetry_state = symmetry_state,
      sig_hull_vol = sig_hull_vol_fmt,
      int_hull_vol = int_hull_vol_fmt,
      vol_ratio = vol_ratio_fmt,
      dist_to_median = dist_to_median_fmt,
      vol_implication = vol_implication
    ) %>%
    arrange(factor(Panel, levels = c("A","B","C","D","E","F")))
  
  out
}

fs1_panel_key <- build_fs1_panel_key(representatives_df)
print(fs1_panel_key, n = Inf)

# ---- Export (TSV is ideal for Supplementary; XLSX optional) ----
out_tsv  <- "Supplementary_Table_FS1_panel_key.tsv"
out_xlsx <- "Supplementary_Table_FS1_panel_key.xlsx"

write.table(fs1_panel_key, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# Optional XLSX export (requires openxlsx)
if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(fs1_panel_key, file = out_xlsx, overwrite = TRUE)
} else {
  message("Package 'openxlsx' not installed; skipping XLSX export. TSV saved: ", out_tsv)
}

################################################################################
# START
################################################################################

###############################################################################
# Figure S1 Final — TRUE 3×2 multiscene HTML (FS1_HTML)
# Uses representatives_df (6 panels, A–F) selected deterministically.
###############################################################################

stopifnot(requireNamespace("plotly", quietly = TRUE))
stopifnot(requireNamespace("htmlwidgets", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("stringr", quietly = TRUE))

library(dplyr)
library(stringr)

# ---- Required objects -------------------------------------------------------
stopifnot(exists("representatives_df"))
stopifnot(all(c("panel","Circuitries_id","vol_implication","complexity_tier","symmetry_state") %in% names(representatives_df)))
stopifnot(exists("tensor"))
stopifnot(exists("embedding"))
stopifnot(exists("build_circuitry_polytope"))
stopifnot(exists("plot_circuitry_polytope_custom"))

# ---------------------------------------------------------------------------
# 0) Output directory
# ---------------------------------------------------------------------------
out_dir_S1 <- "Figure_S1_panels_html"
dir.create(out_dir_S1, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1) Order panels A–F (row-major: Low→Intermediate→High; Symmetric→Sig_dominant)
#    Your representatives_df should already be arranged this way; we enforce it.
# ---------------------------------------------------------------------------
panel_order_S1 <- c("A","B","C","D","E","F")

S1_tagged <- representatives_df %>%
  mutate(panel = as.character(panel)) %>%
  filter(panel %in% panel_order_S1) %>%
  mutate(panel = factor(panel, levels = panel_order_S1)) %>%
  arrange(panel) %>%
  mutate(panel = as.character(panel))

stopifnot(nrow(S1_tagged) == 6)
stopifnot(n_distinct(S1_tagged$panel) == 6)

make_key <- function(x) {
  x %>%
    stringr::str_replace_all("\\s*/\\s*", "_") %>%
    stringr::str_replace_all("\\s+", "_") %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]", "")
}

report_S1 <- S1_tagged %>%
  transmute(
    panel = panel,
    Circuitries_id = as.character(Circuitries_id),
    circuitry_key = make_key(Circuitries_id),
    html_file = file.path(out_dir_S1, paste0("Figure_S1_Panel_", panel, "_", circuitry_key, ".html")),
    complexity_tier = as.character(complexity_tier),
    symmetry_state  = as.character(symmetry_state),
    vol_implication = as.character(vol_implication)
  )

# ---------------------------------------------------------------------------
# 2) Two-line panel titles (collision-safe)
#    Line 1 (bold): Panel X — Circuitries_id
#    Line 2a: complexity_tier | symmetry_state
#    Line 2b: vol_implication
# ---------------------------------------------------------------------------
make_S1_panel_title <- function(panel_label, circuitry_id, complexity_tier, symmetry_state, vol_imp) {
  line1 <- sprintf("Panel %s — %s", panel_label, circuitry_id)
  line2a <- sprintf("%s | %s", complexity_tier, symmetry_state)
  line2b <- sprintf("%s", vol_imp)
  
  paste0(
    "<b>", line1, "</b>",
    "<br><span style='font-size:10px'>", line2a, "</span>",
    "<br><span style='font-size:10px'>", line2b, "</span>"
  )
}

panel_titles_S1 <- vector("list", nrow(report_S1))
for (i in seq_len(nrow(report_S1))) {
  panel_titles_S1[[i]] <- make_S1_panel_title(
    panel_label = report_S1$panel[i],
    circuitry_id = report_S1$Circuitries_id[i],
    complexity_tier = report_S1$complexity_tier[i],
    symmetry_state  = report_S1$symmetry_state[i],
    vol_imp = report_S1$vol_implication[i]
  )
}

# ---------------------------------------------------------------------------
# 3) Latent dimension names (18D) — same as Figure 4 / Figure S2
# ---------------------------------------------------------------------------
latent_dim_names <- c(
  "rho", "rho_strength", "TN_dir", "TN_strength",
  "OS_dir", "OS_strength", "OS_lr_chisq",
  "DSS_dir", "DSS_strength", "DSS_lr_chisq",
  "DFI_dir", "DFI_strength", "DFI_lr_chisq",
  "PFI_dir", "PFI_strength", "PFI_lr_chisq",
  "TME_score", "Immune_dir"
)
stopifnot(length(latent_dim_names) == 18)

# ---------------------------------------------------------------------------
# 4) Build the SIX plotly objects in memory (A–F)
# ---------------------------------------------------------------------------
build_panel_plotly_S1 <- function(circuitry_id, tensor, embedding) {
  
  poly <- build_circuitry_polytope(
    tensor = tensor,
    embedding = embedding,
    circuitry_id = circuitry_id
  )
  
  plot_circuitry_polytope_custom(
    poly_obj = poly,
    title = "",           # no competing title
    outfile_html = NULL
  )
}

fig_list_S1 <- setNames(vector("list", 6), panel_order_S1)
for (k in seq_len(nrow(report_S1))) {
  lbl <- report_S1$panel[k]
  cid <- report_S1$Circuitries_id[k]
  fig_list_S1[[lbl]] <- build_panel_plotly_S1(cid, tensor, embedding)
}

# ---------------------------------------------------------------------------
# 5) TRUE multiscene composer (reused engine; set nrows=3, ncols=2)
# ---------------------------------------------------------------------------
make_multiscene_3d <- function(fig_list,
                               nrows, ncols,
                               margin_x = 0.08, margin_y = 0.08,
                               latent_dim_names = NULL,
                               right_margin_px = 470,
                               top_margin_px = 80,
                               bottom_margin_px = 10,
                               left_margin_px = 10,
                               legend_x = 1.08,
                               legend_y = 0.40,
                               dimblock_x = 1.08,
                               dimblock_y = 0.98,
                               dimblock_fontsize = 11,
                               legend_fontsize = 11,
                               panel_titles = NULL,
                               panel_title_y_offset = 0.004) {
  
  stopifnot(length(fig_list) == nrows * ncols)
  if (!is.null(panel_titles)) stopifnot(length(panel_titles) == length(fig_list))
  
  x_breaks <- seq(0, 1, length.out = ncols + 1)
  y_breaks <- seq(0, 1, length.out = nrows + 1)
  
  x_domains <- lapply(seq_len(ncols), function(j) {
    c(x_breaks[j] + margin_x/2, x_breaks[j+1] - margin_x/2)
  })
  
  y_domains <- lapply(seq_len(nrows), function(i) {
    top_row_index <- nrows - i + 1
    c(y_breaks[top_row_index] + margin_y/2, y_breaks[top_row_index+1] - margin_y/2)
  })
  
  out <- plotly::plot_ly(type = "scatter3d")
  out$x$data <- list()
  
  for (k in seq_along(fig_list)) {
    
    pb <- plotly::plotly_build(fig_list[[k]])
    scene_id <- if (k == 1) "scene" else paste0("scene", k)
    
    for (tr in pb$x$data) {
      tr$scene <- scene_id
      if (k != 1) tr$showlegend <- FALSE
      out$x$data <- c(out$x$data, list(tr))
    }
  }
  
  lay <- list(
    margin = list(l = left_margin_px, r = right_margin_px, b = bottom_margin_px, t = top_margin_px),
    xaxis = list(visible = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(visible = FALSE, showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    legend = list(
      x = legend_x, y = legend_y, xref = "paper", yref = "paper",
      xanchor = "left", yanchor = "top",
      orientation = "v",
      bgcolor = "rgba(255,255,255,0.88)",
      bordercolor = "rgba(0,0,0,0.15)",
      borderwidth = 1,
      font = list(size = legend_fontsize)
    )
  )
  
  for (k in seq_along(fig_list)) {
    
    scene_id <- if (k == 1) "scene" else paste0("scene", k)
    
    row_i <- ceiling(k / ncols)
    col_j <- k - (row_i - 1) * ncols
    
    lay[[scene_id]] <- list(
      domain = list(x = x_domains[[col_j]], y = y_domains[[row_i]]),
      xaxis  = list(title = "PC1"),
      yaxis  = list(title = "PC2"),
      zaxis  = list(title = "PC3")
    )
    
    if (!is.null(panel_titles)) {
      if (is.null(lay$annotations)) lay$annotations <- list()
      lay$annotations <- c(lay$annotations, list(
        list(
          text = panel_titles[[k]],
          x = mean(x_domains[[col_j]]),
          y = y_domains[[row_i]][2] + panel_title_y_offset,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          font = list(size = 12),
          xanchor = "center",
          yanchor = "bottom",
          align = "center"
        )
      ))
    }
  }
  
  if (!is.null(latent_dim_names)) {
    stopifnot(length(latent_dim_names) == 18)
    
    dim_txt <- paste(
      "<b>Latent dimensions (18D):</b>",
      paste(sprintf("%02d: %s", seq_along(latent_dim_names), latent_dim_names), collapse = "<br>"),
      sep = "<br>"
    )
    
    if (is.null(lay$annotations)) lay$annotations <- list()
    lay$annotations <- c(lay$annotations, list(
      list(
        text = dim_txt,
        x = dimblock_x,
        y = dimblock_y,
        xref = "paper",
        yref = "paper",
        xanchor = "left",
        yanchor = "top",
        align = "left",
        showarrow = FALSE,
        bgcolor = "rgba(255,255,255,0.92)",
        bordercolor = "rgba(0,0,0,0.15)",
        borderwidth = 1,
        font = list(size = dimblock_fontsize)
      )
    ))
  }
  
  do.call(plotly::layout, c(list(out), lay))
}

# ---------------------------------------------------------------------------
# 6) Build composite + save + open + export raster/tiff
# ---------------------------------------------------------------------------
figS1_multi <- make_multiscene_3d(
  fig_list = list(fig_list_S1$A, fig_list_S1$B, fig_list_S1$C, fig_list_S1$D, fig_list_S1$E, fig_list_S1$F),
  nrows = 3, ncols = 2,
  
  margin_x = 0.08,
  margin_y = 0.08,
  
  latent_dim_names = latent_dim_names,
  
  right_margin_px = 470,
  legend_x = 1.08,
  dimblock_x = 1.08,
  
  dimblock_y = 0.98,
  legend_y   = 0.40,
  
  top_margin_px = 85,
  panel_titles = panel_titles_S1,
  panel_title_y_offset = 0.006
)

out_file_S1 <- file.path(out_dir_S1, "Figure_S1_MULTISCENE_3x2.html")
htmlwidgets::saveWidget(figS1_multi, out_file_S1, selfcontained = TRUE)

browseURL(paste0("file:///", normalizePath(out_file_S1, winslash = "/")))

# Raster export: increase height because 3 rows
webshot2::webshot(
  out_file_S1,
  file   = "Figure_S1_MULTISCENE_3x2.png",
  vwidth = 2400,
  vheight = 2700,
  zoom   = 2.5
)

magick::image_write(
  magick::image_read("Figure_S1_MULTISCENE_3x2.png"),
  path    = "Figure_S1_MULTISCENE_3x2_600dpi.tiff",
  format  = "tiff",
  density = "600x600"
)

################################################################################
# END
################################################################################

### Creating Dataset_S1 for manuscript
rio::export(circuitries, "Dataset_S1.tsv")

Dataset_S1 <- import("Dataset_S1.tsv")

################################################################################
# START
################################################################################

####
####
####
#### Figure 2 improved TIFF, PNG and PDF
#### Version 1 (V1). Axis scales are real data point coordinates of the 3D polytopes
####
################################################################################
# Figure X – Conceptual Geometry (Panels A–F)
# ------------------------------------------------------------------------------
# LOCKED, REPRODUCIBLE EXPORT PIPELINE (HTML → PNG → 600 dpi TIFF → PDF COMPOSITE)
#
# PURPOSE
#   Generate six synthetic 3D convex-hull panels (A–F) and export publication-
#   grade static rasters while preserving per-panel hull colors and typography.
#
# AUTHORITATIVE DESIGN PRINCIPLES (WHAT THIS SCRIPT ENFORCES)
#   1) Per-panel hull colors are hard-locked (WebGL-safe RGBA) and never depend
#      on default plotly color scales.
#   2) Native 3D axis titles are DISABLED to avoid WebGL label collisions.
#      Axis titles ("Latent axis 1–3") are reintroduced ONLY as PAPER-anchored
#      annotations (collision-proof in static export).
#   3) Bottom tick-label clipping (export artifact) is addressed by allocating
#      paper-space via:
#        - scene$domain (primary lever: lifts/shrinks the 3D scene box), and
#        - margin_b     (secondary lever: increases bottom padding).
#
# IMPORTANT CONSTRAINTS (DO NOT VIOLATE DURING TUNING)
#   - Tick-label font size is fixed by axis_tick and must NOT be changed to
#     resolve clipping. Adjust only domain/margins/annotation placement.
#   - Do NOT introduce camera/range/aspect locking unless explicitly desired,
#     as these can alter perceived geometry (polytope appearance).
#
# CURRENTLY ACTIVE EXPORT/CLIPPING SETTINGS (AS CODED BELOW)
#   - axis titles: disabled in scene; rendered as paper annotations
#   - clipping mitigation: scene_domain_y = c(0.16, 0.92) and margin_b = 40
#   - panel size: 3.5 in × 3.0 in @ 600 dpi (captured via webshot2 viewport)
################################################################################

suppressPackageStartupMessages({
  library(plotly)
  library(htmlwidgets)
  library(geometry)
  library(webshot2)
  library(magick)
})

#-----------------------------
# 1) Build synthetic panels A–F
#-----------------------------
set.seed(42)

generate_point_cloud <- function(type = c("compact","elongated","single_axis","multi_axis","irregular"),
                                 n = 30, scale = 1) {
  type <- match.arg(type)
  
  if (type == "compact") {
    X <- matrix(rnorm(n * 3, sd = 0.35), ncol = 3)
  } else if (type == "elongated") {
    X <- cbind(rnorm(n, sd = 1.20), rnorm(n, sd = 0.25), rnorm(n, sd = 0.25))
  } else if (type == "single_axis") {
    X <- cbind(rnorm(n, sd = 1.50), rnorm(n, sd = 0.12), rnorm(n, sd = 0.12))
  } else if (type == "multi_axis") {
    X <- cbind(rnorm(n, sd = 1.00), rnorm(n, sd = 0.95), rnorm(n, sd = 0.15))
  } else if (type == "irregular") {
    n1 <- floor(n / 3); n2 <- floor(n / 3); n3 <- n - n1 - n2
    X1 <- matrix(rnorm(n1 * 3, mean = -1.2, sd = 0.25), ncol = 3)
    X2 <- matrix(rnorm(n2 * 3, mean =  0.0, sd = 0.30), ncol = 3)
    X3 <- matrix(rnorm(n3 * 3, mean =  1.2, sd = 0.25), ncol = 3)
    X  <- rbind(X1, X2, X3)
  }
  
  X * scale
}

compute_hull <- function(X) {
  ch <- geometry::convhulln(X, options = "FA")
  list(points = X, facets = ch$hull, area = ch$area, vol = ch$vol)
}

plot_hull_3d <- function(hull_obj,
                         title = NULL,
                         point_color = "gray40",
                         point_size  = 8,
                         hull_opacity = 0.35) {
  X <- hull_obj$points
  F <- hull_obj$facets
  stopifnot(is.matrix(X), ncol(X) == 3)
  stopifnot(is.matrix(F), ncol(F) == 3)
  
  plot_ly() %>%
    add_trace(
      x = X[,1], y = X[,2], z = X[,3],
      type = "scatter3d",
      mode = "markers",
      marker = list(size = point_size, color = point_color),
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "mesh3d",
      x = X[,1], y = X[,2], z = X[,3],
      i = F[,1] - 1, j = F[,2] - 1, k = F[,3] - 1,
      opacity = hull_opacity,
      showlegend = FALSE
    ) %>%
    layout(
      title = list(text = title),
      margin = list(l = 0, r = 0, b = 0, t = 20),
      scene = list(
        xaxis = list(title = "Latent axis 1"),
        yaxis = list(title = "Latent axis 2"),
        zaxis = list(title = "Latent axis 3")
      )
    ) %>%
    config(displayModeBar = FALSE)
}

# Build data/hulls
X_A <- generate_point_cloud("compact",     n = 30)
X_B <- generate_point_cloud("compact",     n = 30)
X_C <- generate_point_cloud("elongated",   n = 30)
X_D <- generate_point_cloud("single_axis", n = 30)
X_E <- generate_point_cloud("multi_axis",  n = 30)
X_F <- generate_point_cloud("irregular",   n = 30)

h_A <- compute_hull(X_A); h_B <- compute_hull(X_B); h_C <- compute_hull(X_C)
h_D <- compute_hull(X_D); h_E <- compute_hull(X_E); h_F <- compute_hull(X_F)

fig_A <- plot_hull_3d(h_A, title = "(A) Conceptual projection of a multidimensional omic signature")
fig_B <- plot_hull_3d(h_B, title = "(B) Compact and isotropic geometric organization")
fig_C <- plot_hull_3d(h_C, title = "(C) Elongated and anisotropic geometric organization")
fig_D <- plot_hull_3d(h_D, title = "(D) Coherent single-axis regulatory organization")
fig_E <- plot_hull_3d(h_E, title = "(E) Directional multi-axis regulatory organization")
fig_F <- plot_hull_3d(h_F, title = "(F) Heterogeneous and modular organizational structure")

fig_list <- list(A = fig_A, B = fig_B, C = fig_C, D = fig_D, E = fig_E, F = fig_F)

#-----------------------------
# 2) Per-panel hull colors (LOCKED)
#-----------------------------
hull_palette <- c(
  A = "rgba(243, 156,  18, 0.35)",
  B = "rgba(231,  76,  60, 0.35)",
  C = "rgba(160, 140, 120, 0.35)",
  D = "rgba(140, 140, 140, 0.35)",
  E = "rgba( 52, 152, 219, 0.35)",
  F = "rgba(243, 156,  18, 0.35)"
)

point_color_locked <- "rgba(31, 119, 180, 1.0)"

#-----------------------------
# 3) LOCK STYLE + axis titles as PAPER annotations; 3D axis titles disabled
#    + FIX: prevent bottom tick clipping via scene$domain (robust)
#-----------------------------
lock_panel_style <- function(fig, hull_rgba,
                             point_color  = point_color_locked,
                             point_size   = 10,
                             
                             font_family  = "Arial",
                             base_font    = 24,
                             axis_title   = 50,
                             axis_tick    = 22,
                             
                             # legacy knobs retained (not relied upon)
                             axis_title_standoff = 200,
                             axis_automargin     = TRUE,
                             
                             title_size   = 75,
                             
                             title_y      = 0.85,
                             top_margin_t = 6,
                             
                             # Axis-title annotations (paper coords)
                             xlab_y = 0.02,
                             ylab_x = 0.02,
                             zlab_x = 0.98,
                             
                             # Margins: slightly increased bottom margin to support ticks
                             margin_l = 40,
                             margin_r = 40,
                             margin_b = 40,   # <<< CHANGED (was 80): helps, but domain is the real fix
                             
                             # >>> KEY FIX: lift/shrink the 3D scene so bottom ticks cannot be clipped
                             # Increase domain_y[1] (lower bound) to push scene upward.
                             scene_domain_x = c(0.06, 0.94),
                             scene_domain_y = c(0.16, 0.92),  # <<< CHANGED: lifts scene; prevents bottom tick clipping
                             
                             hull_opacity = 0.35) {
  
  stopifnot(inherits(fig, "plotly"))
  
  pb <- plotly_build(fig)
  if (is.null(pb$x$data) || length(pb$x$data) == 0) {
    stop("plotly_build() produced no traces; invalid plotly object.")
  }
  
  # Extract title text, remove default title
  this_title <- ""
  if (!is.null(pb$x$layout$title) && !is.null(pb$x$layout$title$text)) {
    this_title <- pb$x$layout$title$text
  }
  pb$x$layout$title <- NULL
  
  # Edit traces
  for (i in seq_along(pb$x$data)) {
    tr <- pb$x$data[[i]]
    
    if (identical(tr$type, "scatter3d")) {
      if (is.null(tr$marker)) tr$marker <- list()
      tr$marker$color <- point_color
      tr$marker$size  <- point_size
      pb$x$data[[i]] <- tr
    }
    
    if (identical(tr$type, "mesh3d")) {
      nfaces <- length(tr$i)
      tr$facecolor <- rep(hull_rgba, nfaces)
      tr$opacity   <- hull_opacity
      
      tr$flatshading   <- TRUE
      tr$lighting      <- list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0)
      tr$lightposition <- list(x = 0, y = 0, z = 1e6)
      
      tr$mode  <- NULL
      tr$color <- NULL
      pb$x$data[[i]] <- tr
    }
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # Panel title annotation
  ann_title <- list(
    text      = this_title,
    x         = 0.50,
    y         = title_y,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "bottom",
    showarrow = FALSE,
    font      = list(family = font_family, size = title_size, color = "black")
  )
  
  # Axis-title annotations (paper coords; collision-proof)
  ann_x <- list(
    text      = "Latent axis 1",
    x         = 0.50,
    y         = xlab_y,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "top",
    showarrow = FALSE,
    font      = list(family = font_family, size = axis_title, color = "black")
  )
  
  ann_y <- list(
    text      = "Latent axis 2",
    x         = ylab_x,
    y         = 0.50,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "middle",
    textangle = -90,
    showarrow = FALSE,
    font      = list(family = font_family, size = axis_title, color = "black")
  )
  
  ann_z <- list(
    text      = "Latent axis 3",
    x         = zlab_x,
    y         = 0.50,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "middle",
    textangle = 90,
    showarrow = FALSE,
    font      = list(family = font_family, size = axis_title, color = "black")
  )
  
  pb$x$layout$annotations <- c(pb$x$layout$annotations %||% list(),
                               list(ann_title, ann_x, ann_y, ann_z))
  
  # Disable 3D axis titles; keep tick font unchanged; lift/shrink scene via domain
  pb <- pb %>%
    layout(
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      font = list(family = font_family, size = base_font, color = "black"),
      
      margin = list(l = margin_l, r = margin_r, b = margin_b, t = top_margin_t),
      
      scene = list(
        bgcolor = "white",
        
        # <<< KEY: reserve internal bottom space so tick labels are never clipped
        domain = list(x = scene_domain_x, y = scene_domain_y),
        
        xaxis = list(
          title = list(text = ""),
          tickfont = list(size = axis_tick, family = font_family, color = "black")
        ),
        yaxis = list(
          title = list(text = ""),
          tickfont = list(size = axis_tick, family = font_family, color = "black")
        ),
        zaxis = list(
          title = list(text = ""),
          tickfont = list(size = axis_tick, family = font_family, color = "black")
        )
      )
    ) %>%
    config(displayModeBar = FALSE)
  
  pb
}

# Apply locked style A–F
fig_list_locked <- fig_list
for (nm in names(fig_list_locked)) {
  fig_list_locked[[nm]] <- lock_panel_style(
    fig       = fig_list_locked[[nm]],
    hull_rgba = hull_palette[[nm]]
  )
}

#-----------------------------
# 4) Output folders
#-----------------------------
out_html <- "FigureX_Panels_HTML"
out_png  <- "FigureX_Panels_PNG"
out_tif  <- "FigureX_Panels_TIFF_600dpi"
out_comp <- "FigureX_Composite"
dir.create(out_html, showWarnings = FALSE, recursive = TRUE)
dir.create(out_png,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_tif,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_comp, showWarnings = FALSE, recursive = TRUE)

#-----------------------------
# 5) Export parameters
#-----------------------------
dpi_target <- 600
panel_width_in  <- 3.5
panel_height_in <- 3.0
panel_w_px <- as.integer(round(panel_width_in  * dpi_target))
panel_h_px <- as.integer(round(panel_height_in * dpi_target))

zoom_factor <- 2
render_delay_sec <- 1

#-----------------------------
# 6) Exporter: HTML -> PNG -> TIFF
#-----------------------------
export_plotly_panel <- function(fig, nm,
                                out_html, out_png, out_tif,
                                w_px, h_px,
                                zoom = 2,
                                delay = 1,
                                dpi = 600) {
  
  html_file <- file.path(out_html, sprintf("FigureX_%s.html", nm))
  png_file  <- file.path(out_png,  sprintf("FigureX_%s.png",  nm))
  tif_file  <- file.path(out_tif,  sprintf("FigureX_%s_%ddpi.tif", nm, dpi))
  
  htmlwidgets::saveWidget(fig, file = html_file, selfcontained = FALSE)
  if (!file.exists(html_file)) stop("HTML not created: ", html_file)
  
  webshot2::webshot(
    url     = paste0("file:///", normalizePath(html_file, winslash = "/")),
    file    = png_file,
    vwidth  = w_px,
    vheight = h_px,
    zoom    = zoom,
    delay   = delay
  )
  if (!file.exists(png_file)) stop("PNG not created: ", png_file)
  
  img <- magick::image_read(png_file)
  magick::image_write(
    img,
    path        = tif_file,
    format      = "tiff",
    compression = "lzw",
    density     = paste0(dpi, "x", dpi)
  )
  if (!file.exists(tif_file)) stop("TIFF not created: ", tif_file)
  
  list(html = html_file, png = png_file, tif = tif_file)
}

#-----------------------------
# 7) Export A–F
#-----------------------------
exports <- setNames(vector("list", length(fig_list_locked)), names(fig_list_locked))

for (nm in names(fig_list_locked)) {
  message("\n=== Exporting panel ", nm, " ===")
  exports[[nm]] <- export_plotly_panel(
    fig      = fig_list_locked[[nm]],
    nm       = nm,
    out_html = out_html,
    out_png  = out_png,
    out_tif  = out_tif,
    w_px     = panel_w_px,
    h_px     = panel_h_px,
    zoom     = zoom_factor,
    delay    = render_delay_sec,
    dpi      = dpi_target
  )
  message("OK ", nm, " -> ", exports[[nm]]$tif)
}

#-----------------------------
# 8) Composite 2×3 TIFF + PDF
#-----------------------------
panel_tifs <- vapply(exports, function(x) x$tif, FUN.VALUE = character(1))
panel_tifs <- panel_tifs[c("A","B","C","D","E","F")]

if (any(!file.exists(panel_tifs))) {
  stop("Composite aborted: not all panel TIFF files exist on disk.")
}

imgs <- lapply(panel_tifs, magick::image_read)

pad_px <- 30
imgs <- lapply(imgs, function(im) magick::image_border(im, "white", paste0(pad_px, "x", pad_px)))

row1 <- magick::image_append(magick::image_join(imgs[1:3]), stack = FALSE)
row2 <- magick::image_append(magick::image_join(imgs[4:6]), stack = FALSE)
comp <- magick::image_append(magick::image_join(list(row1, row2)), stack = TRUE)

comp_tif <- file.path(out_comp, sprintf("FigureX_Composite_2x3_%ddpi.tif", dpi_target))
magick::image_write(
  comp,
  path        = comp_tif,
  format      = "tiff",
  compression = "lzw",
  density     = paste0(dpi_target, "x", dpi_target)
)
if (!file.exists(comp_tif)) stop("Composite TIFF not created: ", comp_tif)

comp_pdf_target <- file.path(out_comp, "FigureX_Composite_2x3_V1.pdf")
safe_dir <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
comp_pdf_safe <- file.path(safe_dir, "FigureX_Composite_2x3_V1.pdf")

if (file.exists(comp_pdf_safe)) {
  ok_rm <- tryCatch(file.remove(comp_pdf_safe), error = function(e) FALSE)
  if (!ok_rm && file.exists(comp_pdf_safe)) {
    stop("Cannot overwrite temp PDF (likely open/locked): ", comp_pdf_safe,
         "\nClose the PDF viewer and rerun.")
  }
}

magick::image_write(comp, path = comp_pdf_safe, format = "pdf")
if (!file.exists(comp_pdf_safe)) stop("Composite PDF not created: ", comp_pdf_safe)

ok_copy <- tryCatch(file.copy(comp_pdf_safe, comp_pdf_target, overwrite = TRUE),
                    error = function(e) FALSE)

message("\nDONE.\nPanels written to:\n  HTML: ", normalizePath(out_html),
        "\n  PNG : ", normalizePath(out_png),
        "\n  TIFF: ", normalizePath(out_tif),
        "\nComposite:\n  ", normalizePath(comp_tif),
        if (ok_copy) paste0("\n  ", normalizePath(comp_pdf_target)) else paste0("\n  ", comp_pdf_safe))

###############################################################################
# END
###############################################################################

################################################################################
# START
################################################################################

####
####
####
#### Figure 2 improved TIFF, PNG and PDF
#### Version 2 (V2) axes ticks and numbers are manually placed
####
####
####
################################################################################
# Figure 2 – Conceptual Geometry (Panels A–F)
# ---------------------------------------------------------------------------
# LOCKED, REPRODUCIBLE EXPORT PIPELINE (HTML → PNG → 600 dpi TIFF → PDF COMPOSITE)
# for a SIMULATED / CONCEPTUAL multi-panel 3D geometry figure.

# RATIONALE
#   This figure is illustrative (synthetic point clouds + convex hulls), not
#   quantitative. Therefore, numerical 3D tick labels are not essential to the
#   scientific message, while they are a frequent source of export-time clipping
#   and WebGL rasterization artifacts.

# DESIGN PRINCIPLES (AUTHORITATIVE):
#   1) Native 3D axis TITLES are DISABLED inside the WebGL scene to eliminate
#      label collisions and export-time clipping.
#   2) Axis TITLES ("Latent axis 1–3") are reintroduced exclusively as
#      PAPER-anchored annotations, ensuring collision-free placement.
#   3) Native 3D TICK LABELS are DISABLED entirely (most robust anti-clipping
#      measure for webshot2/WebGL exports).
#   4) A minimal set of tick numbers (conceptual scale cues only) is recreated
#      as PAPER-anchored annotations, avoiding any dependence on WebGL tick
#      rendering.

# IMPORTANT CONSTRAINTS:
#     The 3D geometry itself (camera, scene scaling, point clouds, hulls) must
#     remain unchanged during layout tuning.
#   - Only margins and PAPER annotations are permitted degrees of freedom.

# TARGET OUTPUT:
#     Publication-grade, layout-stable panels suitable for static TIFF/PDF
#     assembly, invariant to browser, device, or WebGL rendering differences.

# NOTE ON SCALE (IMPORTANT):
#   The tick numbers displayed around the panels are *conceptual scale cues only*.
#   They are drawn as PAPER-anchored annotations and are NOT tied to the underlying
#   WebGL axis tick rendering, nor to enforced axis ranges. Therefore, they do NOT
#   necessarily correspond to the true coordinate extents of the simulated 3D point
#   clouds / convex hulls in the Plotly scene.

# NOTE: scene$aspectmode ("data" = data-scaled axes, "cube" = forced equal axes) is intentionally left unset
# (Plotly default) to avoid geometric distortion in this conceptual, simulated figure.

################################################################################

suppressPackageStartupMessages({
  library(plotly)
  library(htmlwidgets)
  library(geometry)
  library(webshot2)
  library(magick)
})

#-----------------------------
# 1) Build synthetic panels A–F
#-----------------------------
set.seed(42)

generate_point_cloud <- function(type = c("compact","elongated","single_axis","multi_axis","irregular"),
                                 n = 30, scale = 1) {
  type <- match.arg(type)
  
  if (type == "compact") {
    X <- matrix(rnorm(n * 3, sd = 0.35), ncol = 3)
  } else if (type == "elongated") {
    X <- cbind(rnorm(n, sd = 1.20), rnorm(n, sd = 0.25), rnorm(n, sd = 0.25))
  } else if (type == "single_axis") {
    X <- cbind(rnorm(n, sd = 1.50), rnorm(n, sd = 0.12), rnorm(n, sd = 0.12))
  } else if (type == "multi_axis") {
    X <- cbind(rnorm(n, sd = 1.00), rnorm(n, sd = 0.95), rnorm(n, sd = 0.15))
  } else if (type == "irregular") {
    n1 <- floor(n / 3); n2 <- floor(n / 3); n3 <- n - n1 - n2
    X1 <- matrix(rnorm(n1 * 3, mean = -1.2, sd = 0.25), ncol = 3)
    X2 <- matrix(rnorm(n2 * 3, mean =  0.0, sd = 0.30), ncol = 3)
    X3 <- matrix(rnorm(n3 * 3, mean =  1.2, sd = 0.25), ncol = 3)
    X  <- rbind(X1, X2, X3)
  }
  
  X * scale
}

compute_hull <- function(X) {
  ch <- geometry::convhulln(X, options = "FA")
  list(points = X, facets = ch$hull, area = ch$area, vol = ch$vol)
}

plot_hull_3d <- function(hull_obj,
                         title = NULL,
                         point_color = "gray40",
                         point_size  = 10,
                         hull_opacity = 0.35) {
  X <- hull_obj$points
  F <- hull_obj$facets
  stopifnot(is.matrix(X), ncol(X) == 3)
  stopifnot(is.matrix(F), ncol(F) == 3)
  
  plot_ly() %>%
    add_trace(
      x = X[,1], y = X[,2], z = X[,3],
      type = "scatter3d",
      mode = "markers",
      marker = list(size = point_size, color = point_color),
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "mesh3d",
      x = X[,1], y = X[,2], z = X[,3],
      i = F[,1] - 1, j = F[,2] - 1, k = F[,3] - 1,
      opacity = hull_opacity,
      showlegend = FALSE
    ) %>%
    layout(
      title = list(text = title),
      margin = list(l = 0, r = 0, b = 0, t = 20),
      scene = list(
        xaxis = list(title = "Latent axis 1"),
        yaxis = list(title = "Latent axis 2"),
        zaxis = list(title = "Latent axis 3")
      )
    ) %>%
    config(displayModeBar = FALSE)
}

# Build data/hulls
X_A <- generate_point_cloud("compact",     n = 30)
X_B <- generate_point_cloud("compact",     n = 30)
X_C <- generate_point_cloud("elongated",   n = 30)
X_D <- generate_point_cloud("single_axis", n = 30)
X_E <- generate_point_cloud("multi_axis",  n = 30)
X_F <- generate_point_cloud("irregular",   n = 30)

h_A <- compute_hull(X_A); h_B <- compute_hull(X_B); h_C <- compute_hull(X_C)
h_D <- compute_hull(X_D); h_E <- compute_hull(X_E); h_F <- compute_hull(X_F)

fig_A <- plot_hull_3d(h_A, title = "(A) Conceptual projection of a multidimensional omic signature")
fig_B <- plot_hull_3d(h_B, title = "(B) Compact and isotropic geometric organization")
fig_C <- plot_hull_3d(h_C, title = "(C) Elongated and anisotropic geometric organization")
fig_D <- plot_hull_3d(h_D, title = "(D) Coherent single-axis regulatory organization")
fig_E <- plot_hull_3d(h_E, title = "(E) Directional multi-axis regulatory organization")
fig_F <- plot_hull_3d(h_F, title = "(F) Heterogeneous and modular organizational structure")

fig_list <- list(A = fig_A, B = fig_B, C = fig_C, D = fig_D, E = fig_E, F = fig_F)

#-----------------------------
# 2) Per-panel hull colors (LOCKED)
#-----------------------------
hull_palette <- c(
  A = "rgba(243, 156,  18, 0.35)",
  B = "rgba(231,  76,  60, 0.35)",
  C = "rgba(160, 140, 120, 0.35)",
  D = "rgba(140, 140, 140, 0.35)",
  E = "rgba( 52, 152, 219, 0.35)",
  F = "rgba(243, 156,  18, 0.35)"
)

point_color_locked <- "rgba(31, 119, 180, 1.0)"

#-----------------------------
# 3) LOCK STYLE + paper axis titles; disable 3D axis titles AND tick labels
#    + recreate conceptual tick numbers as PAPER annotations (export-stable)
#-----------------------------
lock_panel_style <- function(fig, hull_rgba,
                             point_color  = point_color_locked,
                             point_size   = 10,
                             
                             font_family  = "Arial",
                             base_font    = 26,
                             axis_title   = 50,
                             
                             # tick label style (paper-anchored)
                             tick_font    = 50,   # large enough for PDF/TIFF
                             title_size   = 75,
                             title_y      = 0.85,
                             top_margin_t = 20,
                             
                             # Axis-title annotations (paper coords)
                             xlab_y = 0.02,
                             ylab_x = 0.02,
                             zlab_x = 0.98,
                             
                             # Margins
                             margin_l = 40,
                             margin_r = 40,
                             margin_b = 40,
                             
                             # keep your existing domain choice
                             scene_domain_x = c(0.06, 0.94),
                             scene_domain_y = c(0.16, 0.92),
                             
                             hull_opacity = 0.35,
                             
                             # --- PAPER ticks (conceptual scale cues) ---
                             tick_vals = c(-2, 0, 2),
                             
                             # X tick numbers along bottom (paper coords)
                             xtick_x0 = 0.26,
                             xtick_x1 = 0.74,
                             xtick_y  = 0.105,
                             
                             # Y tick numbers along left side (paper coords)
                             ytick_x  = 0.095,
                             ytick_y0 = 0.28,
                             ytick_y1 = 0.72,
                             
                             # Z tick numbers along right side (paper coords)
                             ztick_x  = 0.905,
                             ztick_y0 = 0.28,
                             ztick_y1 = 0.72,
                             
                             # fine shifts to avoid collisions
                             xtick_yshift = 0,
                             ytick_xshift = 0,
                             ztick_xshift = 0) {
  
  stopifnot(inherits(fig, "plotly"))
  
  pb <- plotly_build(fig)
  if (is.null(pb$x$data) || length(pb$x$data) == 0) {
    stop("plotly_build() produced no traces; invalid plotly object.")
  }
  
  # Extract title text, remove default title
  this_title <- ""
  if (!is.null(pb$x$layout$title) && !is.null(pb$x$layout$title$text)) {
    this_title <- pb$x$layout$title$text
  }
  pb$x$layout$title <- NULL
  
  # Edit traces
  for (i in seq_along(pb$x$data)) {
    tr <- pb$x$data[[i]]
    
    if (identical(tr$type, "scatter3d")) {
      if (is.null(tr$marker)) tr$marker <- list()
      tr$marker$color <- point_color
      tr$marker$size  <- point_size
      pb$x$data[[i]] <- tr
    }
    
    if (identical(tr$type, "mesh3d")) {
      nfaces <- length(tr$i)
      tr$facecolor <- rep(hull_rgba, nfaces)
      tr$opacity   <- hull_opacity
      
      tr$flatshading   <- TRUE
      tr$lighting      <- list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0)
      tr$lightposition <- list(x = 0, y = 0, z = 1e6)
      
      tr$mode  <- NULL
      tr$color <- NULL
      pb$x$data[[i]] <- tr
    }
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # Panel title annotation
  ann_title <- list(
    text      = this_title,
    x         = 0.50,
    y         = title_y,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "bottom",
    showarrow = FALSE,
    font      = list(family = font_family, size = title_size, color = "black")
  )
  
  # Axis-title annotations (paper coords)
  ann_x <- list(
    text      = "Latent axis 1",
    x         = 0.50,
    y         = xlab_y,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "top",
    showarrow = FALSE,
    font      = list(family = font_family, size = axis_title, color = "black")
  )
  
  ann_y <- list(
    text      = "Latent axis 2",
    x         = ylab_x,
    y         = 0.50,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "middle",
    textangle = -90,
    showarrow = FALSE,
    font      = list(family = font_family, size = axis_title, color = "black")
  )
  
  ann_z <- list(
    text      = "Latent axis 3",
    x         = zlab_x,
    y         = 0.50,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "middle",
    textangle = 90,
    showarrow = FALSE,
    font      = list(family = font_family, size = axis_title, color = "black")
  )
  
  # -----------------------------
  # PAPER ticks: conceptual ticks at −2, 0, +2 (export-stable)
  # -----------------------------
  # (Scale cues only; not guaranteed to match true 3D coordinates / axis ranges.)
  
  fmt_tick <- function(v) {
    if (abs(v) < 1e-12) return("0")
    if (abs(v - round(v)) < 1e-12) return(as.character(as.integer(round(v))))
    sprintf("%.1f", v)
  }
  
  # >>> AMENDED: positions derived from length(tick_vals) (now 3 ticks)
  x_positions <- seq(xtick_x0, xtick_x1, length.out = length(tick_vals))
  y_positions <- seq(ytick_y0, ytick_y1, length.out = length(tick_vals))
  z_positions <- seq(ztick_y0, ztick_y1, length.out = length(tick_vals))
  
  ann_xticks <- lapply(seq_along(tick_vals), function(k) {
    list(
      text      = fmt_tick(tick_vals[k]),
      x         = x_positions[k],
      y         = xtick_y,
      xref      = "paper",
      yref      = "paper",
      xanchor   = "center",
      yanchor   = "top",
      yshift    = xtick_yshift,
      showarrow = FALSE,
      font      = list(family = font_family, size = tick_font, color = "black")
    )
  })
  
  ann_yticks <- lapply(seq_along(tick_vals), function(k) {
    list(
      text      = fmt_tick(tick_vals[k]),
      x         = ytick_x,
      y         = y_positions[k],
      xref      = "paper",
      yref      = "paper",
      xanchor   = "right",
      yanchor   = "middle",
      xshift    = ytick_xshift,
      showarrow = FALSE,
      font      = list(family = font_family, size = tick_font, color = "black")
    )
  })
  
  ann_zticks <- lapply(seq_along(tick_vals), function(k) {
    list(
      text      = fmt_tick(tick_vals[k]),
      x         = ztick_x,
      y         = z_positions[k],
      xref      = "paper",
      yref      = "paper",
      xanchor   = "left",
      yanchor   = "middle",
      xshift    = ztick_xshift,
      showarrow = FALSE,
      font      = list(family = font_family, size = tick_font, color = "black")
    )
  })
  
  pb$x$layout$annotations <- c(
    list(ann_title, ann_x, ann_y, ann_z),
    ann_xticks,
    ann_yticks,
    ann_zticks
  )
  
  axis3d_no_ticks <- list(
    title = list(text = ""),
    showticklabels = FALSE,
    ticks = "",
    showticksuffix = "none",
    showexponent = "none",
    # --- grid ---
    showgrid       = TRUE,
    gridcolor      = "rgba(0,0,0,0.35)",   # visible gray
    gridwidth      = 2,                    # thicker survives rasterization
    
    # --- axis background planes (critical for grid visibility) ---
    showbackground = TRUE,
    backgroundcolor = "rgba(0,0,0,0.04)",   # very light gray plane (still “white-ish”)
    
    # --- axis line (optional but helps legibility) ---
    showline       = TRUE,
    linecolor      = "rgba(0,0,0,0.35)",
    linewidth      = 2,
    
    # --- zero line ---
    zeroline       = TRUE,
    zerolinecolor  = "rgba(0,0,0,0.45)",
    zerolinewidth  = 2
    # NOTE: range (example: range = c(-2, 2)) intentionally left unset; locking ranges can shrink/deform simulated geometry.
  )
  
  pb <- pb %>%
    layout(
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      font = list(family = font_family, size = base_font, color = "black"),
      
      margin = list(l = margin_l, r = margin_r, b = margin_b, t = top_margin_t),
      
      scene = list(
        bgcolor = "white",
        domain  = list(x = scene_domain_x, y = scene_domain_y),
        # aspectmode ("data" = data-scaled axes, "cube" = forced equal axes) intentionally left unset (Plotly default)
        # to avoid geometric distortion in this conceptual, simulated figure
        xaxis = axis3d_no_ticks,
        yaxis = axis3d_no_ticks,
        zaxis = axis3d_no_ticks
      )
    ) %>%
    config(displayModeBar = FALSE)
  
  pb
}

# Apply locked style A–F
fig_list_locked <- fig_list
for (nm in names(fig_list_locked)) {
  fig_list_locked[[nm]] <- lock_panel_style(
    fig       = fig_list_locked[[nm]],
    hull_rgba = hull_palette[[nm]]
  )
}

#-----------------------------
# 4) Output folders
#-----------------------------
out_html <- "FigureX_Panels_HTML"
out_png  <- "FigureX_Panels_PNG"
out_tif  <- "FigureX_Panels_TIFF_600dpi"
out_comp <- "FigureX_Composite"
dir.create(out_html, showWarnings = FALSE, recursive = TRUE)
dir.create(out_png,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_tif,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_comp, showWarnings = FALSE, recursive = TRUE)

#-----------------------------
# 5) Export parameters
#-----------------------------
dpi_target <- 600
panel_width_in  <- 3.5
panel_height_in <- 3.0
panel_w_px <- as.integer(round(panel_width_in  * dpi_target))
panel_h_px <- as.integer(round(panel_height_in * dpi_target))

zoom_factor <- 2
render_delay_sec <- 1

#-----------------------------
# 6) Exporter: HTML -> PNG -> TIFF
#-----------------------------
export_plotly_panel <- function(fig, nm,
                                out_html, out_png, out_tif,
                                w_px, h_px,
                                zoom = 2,
                                delay = 1,
                                dpi = 600) {
  
  html_file <- file.path(out_html, sprintf("FigureX_%s.html", nm))
  png_file  <- file.path(out_png,  sprintf("FigureX_%s.png",  nm))
  tif_file  <- file.path(out_tif,  sprintf("FigureX_%s_%ddpi.tif", nm, dpi))
  
  htmlwidgets::saveWidget(fig, file = html_file, selfcontained = TRUE)
  if (!file.exists(html_file)) stop("HTML not created: ", html_file)
  
  webshot2::webshot(
    url     = paste0("file:///", normalizePath(html_file, winslash = "/")),
    file    = png_file,
    vwidth  = w_px,
    vheight = h_px,
    zoom    = zoom,
    delay   = delay
  )
  if (!file.exists(png_file)) stop("PNG not created: ", png_file)
  
  img <- magick::image_read(png_file)
  magick::image_write(
    img,
    path        = tif_file,
    format      = "tiff",
    compression = "lzw",
    density     = paste0(dpi, "x", dpi)
  )
  if (!file.exists(tif_file)) stop("TIFF not created: ", tif_file)
  
  list(html = html_file, png = png_file, tif = tif_file)
}

#-----------------------------
# 7) Export A–F
#-----------------------------
exports <- setNames(vector("list", length(fig_list_locked)), names(fig_list_locked))

for (nm in names(fig_list_locked)) {
  message("\n=== Exporting panel ", nm, " ===")
  exports[[nm]] <- export_plotly_panel(
    fig      = fig_list_locked[[nm]],
    nm       = nm,
    out_html = out_html,
    out_png  = out_png,
    out_tif  = out_tif,
    w_px     = panel_w_px,
    h_px     = panel_h_px,
    zoom     = zoom_factor,
    delay    = render_delay_sec,
    dpi      = dpi_target
  )
  message("OK ", nm, " -> ", exports[[nm]]$tif)
}

#-----------------------------
# 8) Composite 2×3 TIFF + PDF
#-----------------------------
panel_tifs <- vapply(exports, function(x) x$tif, FUN.VALUE = character(1))
panel_tifs <- panel_tifs[c("A","B","C","D","E","F")]

if (any(!file.exists(panel_tifs))) {
  stop("Composite aborted: not all panel TIFF files exist on disk.")
}

imgs <- lapply(panel_tifs, magick::image_read)

pad_px <- 30
imgs <- lapply(imgs, function(im) magick::image_border(im, "white", paste0(pad_px, "x", pad_px)))

row1 <- magick::image_append(magick::image_join(imgs[1:3]), stack = FALSE)
row2 <- magick::image_append(magick::image_join(imgs[4:6]), stack = FALSE)
comp <- magick::image_append(magick::image_join(list(row1, row2)), stack = TRUE)

comp_tif <- file.path(out_comp, sprintf("FigureX_Composite_2x3_%ddpi.tif", dpi_target))
magick::image_write(
  comp,
  path        = comp_tif,
  format      = "tiff",
  compression = "lzw",
  density     = paste0(dpi_target, "x", dpi_target)
)
if (!file.exists(comp_tif)) stop("Composite TIFF not created: ", comp_tif)

comp_pdf_target <- file.path(out_comp, "FigureX_Composite_2x3_V2.pdf")
safe_dir <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
comp_pdf_safe <- file.path(safe_dir, "FigureX_Composite_2x3_V2.pdf")

if (file.exists(comp_pdf_safe)) {
  ok_rm <- tryCatch(file.remove(comp_pdf_safe), error = function(e) FALSE)
  if (!ok_rm && file.exists(comp_pdf_safe)) {
    stop("Cannot overwrite temp PDF (likely open/locked): ", comp_pdf_safe,
         "\nClose the PDF viewer and rerun.")
  }
}

magick::image_write(comp, path = comp_pdf_safe, format = "pdf")
if (!file.exists(comp_pdf_safe)) stop("Composite PDF not created: ", comp_pdf_safe)

ok_copy <- tryCatch(file.copy(comp_pdf_safe, comp_pdf_target, overwrite = TRUE),
                    error = function(e) FALSE)

message("\nDONE.\nPanels written to:\n  HTML: ", normalizePath(out_html),
        "\n  PNG : ", normalizePath(out_png),
        "\n  TIFF: ", normalizePath(out_tif),
        "\nComposite:\n  ", normalizePath(comp_tif),
        if (ok_copy) paste0("\n  ", normalizePath(comp_pdf_target)) else paste0("\n  ", comp_pdf_safe))

################################################################################
# END
################################################################################

###############################################################################
# Box 1 | Geometric concepts and definitions used in SigPolytope
# Publication-grade export: PDF (cairo, embedded fonts) + TIFF (ragg, 600 dpi)
# DESIGN UPDATE:
#   (1) Increased inner padding (text farther from box lines)
#   (2) Blue, thicker border
#   (3) Light blue semi-transparent background fill
###############################################################################

suppressPackageStartupMessages({
  library(grid)
})

has_ragg <- requireNamespace("ragg", quietly = TRUE)

box_title <- "Box 1 | Geometric concepts and definitions used in SigPolytope"

box_entries <- list(
  c("Latent space",
    "A shared multidimensional coordinate system derived from the integration of omic, phenotypic, immune, and clinical features. Each latent axis represents a composite dimension capturing correlated biological variation across these layers."),
  c("Signature (sig) and interaction (int) components",
    "For each regulatory circuitry, the sig component represents the upstream regulatory entity (e.g., an ncRNA-based signature), whereas the int component represents the downstream metabolic or molecular interaction target. Each component is represented as a local polytope embedded in latent space."),
  c("Barycenter",
    "The geometric center (mean coordinate) of a polytope, summarizing the net directional contribution of all latent dimensions for a given circuitry component."),
  c("Barycenter distance (d_bary)",
    "The Euclidean distance between the sig and int barycenters in the shared latent space. This metric quantifies the degree of geometric separation between the regulatory and metabolic components of a circuitry."),
  c("Barycenter-distance regimes",
    "Circuitries are discretized into four regimes based on d_bary: high concordance, moderate discordance, strong discordance, and extreme discordance. These regimes provide an interpretable categorization of continuous geometric separation and are used consistently throughout the main text and Supplementary analyses."),
  c("Convex hull",
    "The minimal geometric envelope enclosing all vertices of a polytope. Hull geometry provides a visual and quantitative summary of the internal organization of a signature or interaction component."),
  c("Hull volume",
    "A quantitative measure of the breadth of latent engagement, reflecting the number and diversity of latent dimensions contributing substantially to a given circuitry component."),
  c("Symmetry (sig-int volume ratio)",
    "The relative balance between sig and int hull volumes, indicating whether latent complexity is evenly distributed between components or dominated by one side."),
  c("Separation-with-convergence",
    "A geometric configuration in which sig and int components are spatially separated in latent space (large barycenter distance) yet exhibit concordant sign structure along specific annotated axes (e.g., survival, microenvironmental, immune dimensions). In such cases, barycenter distance captures latent-direction separation, whereas concordance reflects axis-specific functional alignment.")
)

# ---- Page settings ------------------------------------------------------------
page_w_in <- 7
page_h_in <- 9

# ---- Box styling (NEW) --------------------------------------------------------
border_col <- "#1F5AA6"       # deep publication blue
border_lwd <- 2.0            # thicker lines

# Light blue background with transparency (alpha):
# use rgb() so it works consistently across devices
fill_col <- rgb(0.85, 0.92, 1.00, alpha = 0.55)

# ---- Layout / padding (NEW) ---------------------------------------------------
# Outer box margins in npc (position of the box on the page)
outer_left   <- 0.06
outer_right  <- 0.94
outer_bottom <- 0.06
outer_top    <- 0.94

# Inner padding: the text area starts further inside the box
# Increase these to push text farther from the border.
pad_x <- 0.035   # horizontal padding inside box
pad_y <- 0.035   # vertical padding inside box

# Derived text region (inside the box with padding)
text_left   <- outer_left + pad_x
text_right  <- outer_right - pad_x
text_bottom <- outer_bottom + pad_y
text_top    <- outer_top - pad_y

# Typography
title_fs <- 13
body_fs  <- 10.5

# Text wrapping width (characters per line; tune if needed)
wrap_width <- 88

# Vertical spacing controls
lineheight <- 1.15
gap_after_title_nlines <- 0.9
gap_between_entries_nlines <- 0.55

draw_box1 <- function(family = "Helvetica") {
  
  grid.newpage()
  
  # Outer box (blue border + light blue semi-transparent fill)
  grid.rect(
    x = outer_left, y = outer_bottom,
    width  = (outer_right - outer_left),
    height = (outer_top - outer_bottom),
    just = c("left", "bottom"),
    gp = gpar(col = border_col, fill = fill_col, lwd = border_lwd)
  )
  
  # Convert line units to npc height
  points_per_in <- 72
  body_line_pt <- body_fs * lineheight
  page_h_pt <- page_h_in * points_per_in
  line_npc <- body_line_pt / page_h_pt
  
  # Title (now positioned using padded text region, not box edge)
  grid.text(
    box_title,
    x = text_left, y = text_top,
    just = c("left", "top"),
    gp = gpar(fontfamily = family, fontface = "bold", fontsize = title_fs)
  )
  
  # Start y for body text under title (still within padded region)
  y <- text_top - ((title_fs * 1.3) / page_h_pt) - gap_after_title_nlines * line_npc
  
  for (entry in box_entries) {
    term <- entry[1]
    defn <- entry[2]
    
    # Term (bold)
    grid.text(
      term,
      x = text_left, y = y,
      just = c("left", "top"),
      gp = gpar(fontfamily = family, fontface = "bold", fontsize = body_fs)
    )
    y <- y - 1.05 * line_npc
    
    # Wrapped definition lines
    def_lines <- strwrap(defn, width = wrap_width)
    
    for (ln in def_lines) {
      grid.text(
        ln,
        x = text_left, y = y,
        just = c("left", "top"),
        gp = gpar(fontfamily = family, fontface = "plain", fontsize = body_fs)
      )
      y <- y - 1.0 * line_npc
    }
    
    # Gap between entries
    y <- y - gap_between_entries_nlines * line_npc
    
    # Safety: stop if we run out of vertical space in padded region
    if (y < text_bottom + 0.01) {
      warning("Box content exceeded padded page height. Reduce body_fs/wrap_width or increase page_h_in.")
      break
    }
  }
}

# ---- Export -------------------------------------------------------------------
# PDF (vector, embedded fonts)
grDevices::cairo_pdf("Box1_SigPolytope.pdf", width = page_w_in, height = page_h_in, family = "Helvetica")
draw_box1(family = "Helvetica")
dev.off()

# TIFF (600 dpi)
if (has_ragg) {
  ragg::agg_tiff(
    "Box1_SigPolytope_600dpi.tiff",
    width = page_w_in, height = page_h_in, units = "in",
    res = 600, compression = "lzw"
  )
  draw_box1(family = "Helvetica")
  dev.off()
} else {
  tiff(
    "Box1_SigPolytope_600dpi.tiff",
    width = page_w_in, height = page_h_in, units = "in",
    res = 600, compression = "lzw"
  )
  draw_box1(family = "Helvetica")
  dev.off()
}

message("Exported: Box1_SigPolytope.pdf and Box1_SigPolytope_600dpi.tiff")


###############################################################################
# Box 1 | Geometric concepts and definitions used in SigPolytope
# Publication-ready export as:
#   (1) PDF (vector)  : Box1_SigPolytope.pdf
#   (2) TIFF (600 dpi): Box1_SigPolytope_600dpi.tiff
#
# Design requirements implemented:
#   (1) generous padding: text away from ALL lines
#   (2) thick BLUE border
#   (3) light-blue translucent background in body
#   (4) DARKER blue title banner with WHITE title text
#
# Robustness:
#   - avoids gridtext/richtext_grob (previous API mismatch issues)
#   - uses base grid + strwrap for predictable wrapping
#   - TIFF uses ragg if available; otherwise falls back to base tiff()
###############################################################################

suppressPackageStartupMessages({
  library(grid)
})

# ---- 1) Box content ----------------------------------------------------------

box_title <- "Box 1 | Geometric concepts and definitions used in SigPolytope"

box_entries <- list(
  list(
    term = "Latent space",
    def  = "A shared multidimensional coordinate system derived from the integration of omic, phenotypic, immune, and clinical features. Each latent axis represents a composite dimension capturing correlated biological variation across these layers."
  ),
  list(
    term = "Signature (sig) and interaction (int) components",
    def  = "For each regulatory circuitry, the sig component represents the upstream regulatory entity (e.g., an ncRNA-based signature), whereas the int component represents the downstream metabolic or molecular interaction target. Each component is represented as a local polytope embedded in latent space."
  ),
  list(
    term = "Barycenter",
    def  = "The geometric center (mean coordinate) of a polytope, summarizing the net directional contribution of all latent dimensions for a given circuitry component."
  ),
  list(
    term = "Barycenter distance (d_bary)",
    def  = "The Euclidean distance between the sig and int barycenters in the shared latent space. This metric quantifies the degree of geometric separation between the regulatory and metabolic components of a circuitry."
  ),
  list(
    term = "Barycenter-distance regimes",
    def  = "Circuitries are discretized into four regimes based on d_bary—high concordance, moderate discordance, strong discordance, and extreme discordance. These regimes provide an interpretable categorization of continuous geometric separation and are used consistently throughout the main text and Supplementary analyses."
  ),
  list(
    term = "Convex hull",
    def  = "The minimal geometric envelope enclosing all vertices of a polytope. Hull geometry provides a visual and quantitative summary of the internal organization of a signature or interaction component."
  ),
  list(
    term = "Hull volume",
    def  = "A quantitative measure of the breadth of latent engagement, reflecting the number and diversity of latent dimensions contributing substantially to a given circuitry component."
  ),
  list(
    term = "Symmetry (sig–int volume ratio)",
    def  = "The relative balance between sig and int hull volumes, indicating whether latent complexity is evenly distributed between components or dominated by one side."
  ),
  list(
    term = "Separation-with-convergence",
    def  = "A geometric configuration in which sig and int components are spatially separated in latent space (large barycenter distance) yet exhibit concordant sign structure along specific annotated axes (e.g., survival, microenvironmental, or immune dimensions). In such cases, barycenter distance captures latent-direction separation, whereas concordance reflects axis-specific functional alignment."
  )
)

# ---- 2) Visual design parameters --------------------------------------------

# Page size (inches) – adjust if your journal has strict box sizing
page_w_in <- 7
page_h_in <- 9

# Fonts (Helvetica is safe and avoids “font weirdness”)
family <- "Helvetica"

# Border + fills
border_col   <- "#1F4FB2"                 # blue border
border_lwd   <- 2.2                       # thick
body_fill    <- grDevices::adjustcolor("#D9E8FF", alpha.f = 0.45)  # light blue transparent
title_bg_col <- "#1F4FB2"                 # darker blue title banner
title_col    <- "white"

# Layout fractions (NPC coordinates)
outer_left   <- 0.06
outer_right  <- 0.94
outer_bottom <- 0.06
outer_top    <- 0.94

title_height_frac <- 0.10                 # title banner height fraction (of box height)

# Padding inside body area (NPC)
pad_x <- 0.03
pad_y <- 0.035

# Typography
title_fs <- 13
term_fs  <- 10.6
def_fs   <- 10.2

# Line spacing controls (NPC)
lineheight_term <- 1.15
lineheight_def  <- 1.25

# Spacing between term and definition (in “lines”)
gap_term_to_def_lines <- 0.25
gap_between_entries_lines <- 0.60

# ---- 3) Helper: wrap text to a character width --------------------------------
# We use a conservative approximation to keep wrapping stable across devices.
approx_wrap_width_chars <- function(body_width_npc, fontsize_pt, page_w_in) {
  # Approximate characters per line:
  # typical glyph width ~ 0.52 * fontsize (pt) in points
  # points per inch = 72
  body_width_in <- body_width_npc * page_w_in
  body_width_pt <- body_width_in * 72
  chars <- floor(body_width_pt / (0.52 * fontsize_pt))
  max(55, min(chars, 120))  # clamp to safe range
}

# ---- 4) Renderer -------------------------------------------------------------

draw_box1 <- function() {
  
  grid.newpage()
  
  box_w <- outer_right - outer_left
  box_h <- outer_top - outer_bottom
  
  # --- Body rectangle (light blue translucent) ---
  grid.rect(
    x = outer_left,
    y = outer_bottom,
    width  = box_w,
    height = box_h * (1 - title_height_frac),
    just = c("left", "bottom"),
    gp = gpar(col = border_col, fill = body_fill, lwd = border_lwd)
  )
  
  # --- Title banner rectangle (dark blue) ---
  grid.rect(
    x = outer_left,
    y = outer_bottom + box_h * (1 - title_height_frac),
    width  = box_w,
    height = box_h * title_height_frac,
    just = c("left", "bottom"),
    gp = gpar(col = border_col, fill = title_bg_col, lwd = border_lwd)
  )
  
  # --- Title text (white, centered vertically in banner) ---
  grid.text(
    box_title,
    x = outer_left + pad_x,
    y = outer_bottom + box_h * (1 - title_height_frac/2),
    just = c("left", "center"),
    gp = gpar(fontfamily = family, fontface = "bold", fontsize = title_fs, col = title_col)
  )
  
  # --- Body text region (viewport) ---
  body_left   <- outer_left + pad_x
  body_right  <- outer_right - pad_x
  body_top    <- outer_bottom + box_h * (1 - title_height_frac) - pad_y
  body_bottom <- outer_bottom + pad_y
  
  body_w_npc <- body_right - body_left
  wrap_term  <- approx_wrap_width_chars(body_w_npc, term_fs, page_w_in)
  wrap_def   <- approx_wrap_width_chars(body_w_npc, def_fs,  page_w_in)
  
  # Convert “one line” height in NPC (approx) from point size
  line_npc_term <- (term_fs / (page_h_in * 72)) * lineheight_term
  line_npc_def  <- (def_fs  / (page_h_in * 72)) * lineheight_def
  
  # Start cursor at top of body
  y <- body_top
  
  # Draw each entry
  for (e in box_entries) {
    
    # Term (bold), wrapped
    term_lines <- strwrap(e$term, width = wrap_term)
    for (ln in term_lines) {
      grid.text(
        ln,
        x = body_left, y = y,
        just = c("left", "top"),
        gp = gpar(fontfamily = family, fontface = "bold", fontsize = term_fs, col = "black")
      )
      y <- y - line_npc_term
      if (y < body_bottom) break
    }
    
    # Small gap term -> definition
    y <- y - gap_term_to_def_lines * line_npc_def
    if (y < body_bottom) break
    
    # Definition (regular), wrapped
    def_lines <- strwrap(e$def, width = wrap_def)
    for (ln in def_lines) {
      grid.text(
        ln,
        x = body_left, y = y,
        just = c("left", "top"),
        gp = gpar(fontfamily = family, fontface = "plain", fontsize = def_fs, col = "black")
      )
      y <- y - line_npc_def
      if (y < body_bottom) break
    }
    
    # Gap between entries
    y <- y - gap_between_entries_lines * line_npc_def
    if (y < body_bottom) break
  }
  
  invisible(TRUE)
}

# ---- 5) Export helpers -------------------------------------------------------

export_box1_pdf <- function(filename = "Box1_SigPolytope.pdf") {
  grDevices::pdf(
    file = filename,
    width = page_w_in,
    height = page_h_in,
    useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_box1()
  message("PDF saved: ", normalizePath(filename))
}

export_box1_tiff <- function(filename = "Box1_SigPolytope_600dpi.tiff", dpi = 600) {
  # Prefer ragg for high-quality text rasterization
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(
      filename = filename,
      width  = page_w_in,
      height = page_h_in,
      units  = "in",
      res    = dpi,
      compression = "lzw"
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_box1()
    message("TIFF saved (ragg): ", normalizePath(filename))
  } else {
    # Fallback: base tiff device
    grDevices::tiff(
      filename = filename,
      width  = page_w_in,
      height = page_h_in,
      units  = "in",
      res    = dpi,
      compression = "lzw",
      type = "cairo"
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_box1()
    message("TIFF saved (base tiff): ", normalizePath(filename))
  }
}

# ---- 6) Run exports ----------------------------------------------------------

export_box1_pdf("Box1_SigPolytope.pdf")
export_box1_tiff("Box1_SigPolytope_600dpi.tiff", dpi = 600)

##### 
##### 
##### Box 1. Rounded lines version, italics and subscript
##### 
##### 
###############################################################################
# Box 1 | Geometric concepts and definitions used in SigPolytope
# Publication-ready export as PDF (vector) and TIFF (600 dpi)
#
# This revision CERTIFIES formatting globally:
#  - Every "d_bary" or "dbary" in ALL descriptions renders as italic(d)[bary]
#  - Every "sig" and "int" in ALL descriptions renders in italics
#
# Implementation:
#  - Terms remain plotmath expressions (bold + italic/subscript where needed)
#  - Descriptions remain wrapped + justified EXCEPT lines containing tokens
#    {d_bary, dbary, sig, int}, which are token-typeset (no justification)
#  - Box is centered horizontally + vertically on page (equal margins)
#  - Auto fit-to-box prevents bottom clipping
###############################################################################

suppressPackageStartupMessages({
  library(grid)
})

have_ragg  <- requireNamespace("ragg", quietly = TRUE)
have_cairo <- capabilities("cairo")

# ---- 0) Text sanitization -----------------------------------------------------
sanitize_ascii <- function(x) {
  x <- gsub("\u2013|\u2014", "-", x, perl = TRUE)
  x <- gsub("\u2212", "-", x, perl = TRUE)
  x <- gsub("\u00A0", " ", x, perl = TRUE)
  iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
}

# ---- 1) Content ---------------------------------------------------------------
box_title <- sanitize_ascii("Box 1 | Geometric concepts and definitions used in SigPolytope")

box_entries <- list(
  list(
    term = expression(bold("Latent space")),
    desc = "A shared multidimensional coordinate system derived from the integration of omic, phenotypic, immune, and clinical features. Each latent axis represents a composite dimension capturing correlated biological variation across these layers."
  ),
  
  # ---- NEW (semantic disambiguation anchor; no styling changes) ---------------
  list(
    term = expression(bold("Principal axes")),
    desc = "The orthogonal dimensions of the latent embedding that capture dominant sources of variance across integrated omic, phenotypic, immune, and clinical features."
  ),
  list(
    term = expression(bold("Anisotropy")),
    desc = "A geometric property of the convex hull describing unequal extension of a signature or circuitry across latent axes, reflecting directional dominance rather than positional variance."
  ),
  # ---------------------------------------------------------------------------
  
  list(
    term = expression(bold("Signature ("*italic(sig)*") and interaction ("*italic(int)*") components")),
    desc = "For each regulatory circuitry, the sig component represents the upstream regulatory entity (e.g., an ncRNA-based signature), whereas the int component represents the downstream metabolic or molecular interaction target. Each component is represented as a local polytope embedded in latent space."
  ),
  list(
    term = expression(bold("Barycenter")),
    desc = "The geometric center (mean coordinate) of a polytope, summarizing the net directional contribution of all latent dimensions for a given circuitry component."
  ),
  list(
    term = expression(bold("Barycenter distance ("*italic(d)[bary]*")")),
    desc = "The Euclidean distance between the sig and int barycenters in the shared latent space. This metric quantifies the degree of geometric separation between the regulatory and metabolic components of a circuitry."
  ),
  list(
    term = expression(bold("Barycenter-distance regimes")),
    desc = "Circuitries are discretized into four regimes based on d_bary - high concordance, moderate discordance, strong discordance, and extreme discordance. These regimes provide an interpretable categorization of continuous geometric separation and are used consistently throughout the main text and Supplementary analyses."
  ),
  list(
    term = expression(bold("Convex hull")),
    desc = "The minimal geometric envelope enclosing all vertices of a polytope. Hull geometry provides a visual and quantitative summary of the internal organization of a signature or interaction component."
  ),
  list(
    term = expression(bold("Hull volume")),
    desc = "A quantitative measure of the breadth of latent engagement, reflecting the number and diversity of latent dimensions contributing substantially to a given circuitry component."
  ),
  list(
    term = expression(bold("Symmetry ("*italic(sig)*"-"*italic(int)*" volume ratio)")),
    desc = "The relative balance between sig and int hull volumes, indicating whether latent complexity is evenly distributed between components or dominated by one side."
  ),
  list(
    term = expression(bold("Separation-with-convergence")),
    desc = "A geometric configuration in which sig and int components are spatially separated in latent space (large barycenter distance) yet exhibit concordant sign structure along specific annotated axes (e.g., survival, microenvironmental, immune). In such cases, barycenter distance captures latent-direction separation, whereas concordance reflects axis-specific functional alignment."
  )
)


for (i in seq_along(box_entries)) box_entries[[i]]$desc <- sanitize_ascii(box_entries[[i]]$desc)

# ---- 2) Layout controls -------------------------------------------------------
page_w_in <- 7
page_h_in <- 9
font_family <- "Helvetica"

border_col <- "#1F4E79"
title_fill_alpha  <- adjustcolor("#2F6FAE", alpha.f = 0.35)
title_fill_opaque <- "#2F6FAE"
body_fill          <- adjustcolor("#DCEBFA", alpha.f = 0.55)

corner_r <- unit(7, "pt")
lwd_box  <- 2.4

# centered composite box geometry (npc)
w       <- 0.88
h_body  <- 0.76
h_title <- 0.055
total_h <- h_body + h_title

x0 <- (1 - w) / 2
y_body0  <- (1 - total_h) / 2
y_title0 <- y_body0 + h_body

# padding inside body
pad_x <- 0.04
pad_y <- 0.028

# base typography (auto-fit)
fs_title <- 13.5
fs_term_base <- 10.8
fs_desc_base <- 10.5

term_lh_mult <- 1.12
desc_lh_mult <- 1.10
gap_mult     <- 0.55

target_bottom_gap_pt <- 5

# ---- 3) Wrapping --------------------------------------------------------------
wrap_to_unit_width <- function(text, max_width, gp) {
  text <- gsub("\\s+", " ", trimws(text))
  if (!nzchar(text)) return(character(0))
  
  words <- strsplit(text, " ", fixed = TRUE)[[1]]
  lines <- character(0)
  cur <- ""
  
  max_w_in <- convertWidth(max_width, "in", valueOnly = TRUE)
  for (wd in words) {
    cand <- if (cur == "") wd else paste(cur, wd)
    cand_w <- convertWidth(grobWidth(textGrob(cand, gp = gp)), "in", valueOnly = TRUE)
    if (cand_w <= max_w_in) {
      cur <- cand
    } else {
      if (nzchar(cur)) lines <- c(lines, cur)
      cur <- wd
    }
  }
  if (nzchar(cur)) lines <- c(lines, cur)
  lines
}

# ---- 4) Justification helper --------------------------------------------------
justify_line_to_width <- function(line, target_width, gp) {
  line <- gsub("\\s+", " ", trimws(line))
  if (!nzchar(line)) return(line)
  
  words <- strsplit(line, " ", fixed = TRUE)[[1]]
  if (length(words) < 2) return(line)
  
  cur_w_in <- convertWidth(grobWidth(textGrob(line, gp = gp)), "in", valueOnly = TRUE)
  tgt_w_in <- convertWidth(target_width, "in", valueOnly = TRUE)
  if (cur_w_in >= tgt_w_in) return(line)
  
  sp_w_in <- convertWidth(grobWidth(textGrob(" ", gp = gp)), "in", valueOnly = TRUE)
  if (!is.finite(sp_w_in) || sp_w_in <= 0) return(line)
  
  extra_in <- tgt_w_in - cur_w_in
  extra_spaces <- floor(extra_in / sp_w_in)
  if (extra_spaces <= 0) return(line)
  
  gaps <- length(words) - 1
  add_each <- extra_spaces %/% gaps
  rem      <- extra_spaces %% gaps
  
  out <- words[1]
  for (i in seq_len(gaps)) {
    add_i <- add_each + if (i <= rem) 1 else 0
    out <- paste0(out, paste(rep(" ", 1 + add_i), collapse = ""), words[i + 1])
  }
  out
}

# ---- 5) Token typesetting for descriptions -----------------------------------
# Enforces:
#  - sig, int -> italic
#  - d_bary, dbary -> italic(d)[bary]
needs_token_typeset <- function(line) {
  grepl("\\b(sig|int)\\b|\\b(d_bary|dbary)\\b", line, perl = TRUE)
}

# Draw a single line left-to-right, preserving punctuation and spaces.
draw_line_with_tokens <- function(line, x, y, gp_plain, gp_italic, gp_math) {
  
  m <- gregexpr("\\s+|\\S+", line, perl = TRUE)[[1]]
  tokens <- regmatches(line, list(m))[[1]]
  
  cur_x <- x
  
  for (tok in tokens) {
    
    # spaces
    if (grepl("^\\s+$", tok)) {
      cur_x <- cur_x + grobWidth(textGrob(tok, gp = gp_plain))
      next
    }
    
    # Separate leading/trailing punctuation from the core token
    # e.g. "sig," -> lead="" core="sig" trail=","
    lead  <- sub("^([\\(\\[\\{\\\"\\'\\<\\>\\.,;:!?-]*).*$", "\\1", tok, perl = TRUE)
    trail <- sub("^.*?([\\)\\]\\}\\\"\\'\\<\\>\\.,;:!?-]*)$", "\\1", tok, perl = TRUE)
    core  <- tok
    core  <- sub("^([\\(\\[\\{\\\"\\'\\<\\>\\.,;:!?-]*)", "", core, perl = TRUE)
    core  <- sub("([\\)\\]\\}\\\"\\'\\<\\>\\.,;:!?-]*)$", "", core, perl = TRUE)
    if (!nzchar(core)) core <- tok  # fallback safety
    
    # draw leading punctuation (plain)
    if (nzchar(lead)) {
      grid.draw(textGrob(lead, x = cur_x, y = y, just = c("left", "top"), gp = gp_plain))
      cur_x <- cur_x + grobWidth(textGrob(lead, gp = gp_plain))
    }
    
    # draw core token with enforced formatting
    if (core %in% c("sig", "int")) {
      grid.draw(textGrob(core, x = cur_x, y = y, just = c("left", "top"), gp = gp_italic))
      cur_x <- cur_x + grobWidth(textGrob(core, gp = gp_italic))
      
    } else if (core %in% c("d_bary", "dbary")) {
      # render as italic(d)[bary]
      grid.draw(textGrob(expression(italic(d)[bary]),
                         x = cur_x, y = y, just = c("left", "top"), gp = gp_math))
      cur_x <- cur_x + grobWidth(textGrob(expression(italic(d)[bary]), gp = gp_math))
      
    } else {
      grid.draw(textGrob(core, x = cur_x, y = y, just = c("left", "top"), gp = gp_plain))
      cur_x <- cur_x + grobWidth(textGrob(core, gp = gp_plain))
    }
    
    # draw trailing punctuation (plain)
    if (nzchar(trail)) {
      grid.draw(textGrob(trail, x = cur_x, y = y, just = c("left", "top"), gp = gp_plain))
      cur_x <- cur_x + grobWidth(textGrob(trail, gp = gp_plain))
    }
  }
}

# ---- 6) Build layout ----------------------------------------------------------
build_body_layout <- function(avail_w, gp_desc, lh_term_pt, lh_desc_pt, gap_pt) {
  layout <- list()
  n_terms <- 0; n_desc <- 0; n_gaps <- 0
  
  for (i in seq_along(box_entries)) {
    
    e <- box_entries[[i]]
    layout[[length(layout) + 1]] <- list(type = "term", label = e$term)
    n_terms <- n_terms + 1
    
    desc_lines <- wrap_to_unit_width(e$desc, avail_w, gp_desc)
    if (length(desc_lines) >= 1) {
      for (j in seq_along(desc_lines)) {
        is_last <- (j == length(desc_lines))
        line <- desc_lines[j]
        
        layout[[length(layout) + 1]] <- list(
          type = "desc",
          text = line,
          # do NOT justify token-typeset lines (cannot be safely space-expanded)
          justify = (!is_last) && !needs_token_typeset(line),
          token_typeset = needs_token_typeset(line)
        )
        n_desc <- n_desc + 1
      }
    }
    
    if (i < length(box_entries)) {
      layout[[length(layout) + 1]] <- list(type = "gap")
      n_gaps <- n_gaps + 1
    }
  }
  
  list(
    layout = layout,
    total_height_pt = (n_terms * lh_term_pt) + (n_desc * lh_desc_pt) + (n_gaps * gap_pt)
  )
}

# ---- 7) Renderer --------------------------------------------------------------
draw_box1 <- function(title_fill) {
  
  grid.newpage()
  
  # Title box
  grid.roundrect(
    x = x0, y = y_title0, width = w, height = h_title,
    just = c("left", "bottom"),
    r = corner_r,
    gp = gpar(col = border_col, fill = title_fill, lwd = lwd_box)
  )
  
  # Body box
  grid.roundrect(
    x = x0, y = y_body0, width = w, height = h_body,
    just = c("left", "bottom"),
    r = corner_r,
    gp = gpar(col = border_col, fill = body_fill, lwd = lwd_box)
  )
  
  # Title text centered
  pushViewport(viewport(x = x0, y = y_title0, width = w, height = h_title,
                        just = c("left", "bottom"), clip = "on"))
  grid.text(
    box_title,
    x = unit(0.5, "npc"), y = unit(0.5, "npc"),
    just = "center",
    gp = gpar(fontsize = fs_title, fontface = "bold", col = "black", fontfamily = font_family)
  )
  popViewport()
  
  # Body viewport
  pushViewport(viewport(
    x = x0 + pad_x, y = y_body0 + pad_y,
    width  = w - 2 * pad_x,
    height = h_body - 2 * pad_y,
    just = c("left", "bottom"),
    clip = "on"
  ))
  
  avail_w <- unit(1, "npc")
  avail_h_pt <- convertHeight(unit(1, "npc"), "pt", valueOnly = TRUE)
  
  # Auto fit-to-box scaling
  scale <- 1.0
  scale_min <- 0.90
  
  repeat {
    fs_term <- fs_term_base * scale
    fs_desc <- fs_desc_base * scale
    lh_term_pt <- fs_term * term_lh_mult
    lh_desc_pt <- fs_desc * desc_lh_mult
    gap_pt     <- fs_desc * gap_mult
    
    gp_desc_plain <- gpar(fontface = "plain", fontsize = fs_desc, col = "black", fontfamily = font_family)
    layout_obj <- build_body_layout(avail_w, gp_desc_plain, lh_term_pt, lh_desc_pt, gap_pt)
    
    if (layout_obj$total_height_pt <= (avail_h_pt - target_bottom_gap_pt)) break
    if (scale <= scale_min) break
    scale <- max(scale_min, scale * 0.985)
  }
  
  # Final GPs
  fs_term <- fs_term_base * scale
  fs_desc <- fs_desc_base * scale
  lh_term_pt <- fs_term * term_lh_mult
  lh_desc_pt <- fs_desc * desc_lh_mult
  gap_pt     <- fs_desc * gap_mult
  
  gp_term <- gpar(fontsize = fs_term, col = "black", fontfamily = font_family)
  gp_desc_plain  <- gpar(fontface = "plain",  fontsize = fs_desc, col = "black", fontfamily = font_family)
  gp_desc_italic <- gpar(fontface = "italic", fontsize = fs_desc, col = "black", fontfamily = font_family)
  # gp_math: do not force italic here; expression already includes italic(d)
  gp_math <- gpar(fontface = "plain", fontsize = fs_desc, col = "black", fontfamily = font_family)
  
  layout_obj <- build_body_layout(avail_w, gp_desc_plain, lh_term_pt, lh_desc_pt, gap_pt)
  
  top_offset_pt <- max(0, avail_h_pt - layout_obj$total_height_pt - target_bottom_gap_pt)
  
  cur_x <- unit(0, "npc")
  cur_y <- unit(1, "npc") - unit(top_offset_pt, "pt")
  
  for (item in layout_obj$layout) {
    
    if (identical(item$type, "term")) {
      grid.draw(textGrob(item$label, x = cur_x, y = cur_y, just = c("left", "top"), gp = gp_term))
      cur_y <- cur_y - unit(lh_term_pt, "pt")
      
    } else if (identical(item$type, "desc")) {
      
      if (isTRUE(item$token_typeset)) {
        draw_line_with_tokens(item$text, x = cur_x, y = cur_y,
                              gp_plain = gp_desc_plain, gp_italic = gp_desc_italic, gp_math = gp_math)
      } else {
        ln <- item$text
        if (isTRUE(item$justify)) ln <- justify_line_to_width(ln, target_width = avail_w, gp = gp_desc_plain)
        grid.text(ln, x = cur_x, y = cur_y, just = c("left", "top"), gp = gp_desc_plain)
      }
      
      cur_y <- cur_y - unit(lh_desc_pt, "pt")
      
    } else if (identical(item$type, "gap")) {
      cur_y <- cur_y - unit(gap_pt, "pt")
    }
  }
  
  popViewport()
}

# ---- 8) Export ----------------------------------------------------------------
pdf_file  <- "Box1_SigPolytope_rounded_lines.pdf"
tiff_file <- "Box1_SigPolytope_600dpi_rounded_lines.tiff"

if (have_cairo) {
  grDevices::cairo_pdf(filename = pdf_file, width = page_w_in, height = page_h_in, family = font_family)
  draw_box1(title_fill = title_fill_alpha)
  dev.off()
} else {
  grDevices::pdf(file = pdf_file, width = page_w_in, height = page_h_in,
                 family = font_family, useDingbats = FALSE)
  draw_box1(title_fill = title_fill_opaque)
  dev.off()
}

if (have_ragg) {
  ragg::agg_tiff(filename = tiff_file, width = page_w_in, height = page_h_in,
                 units = "in", res = 600, compression = "lzw")
  draw_box1(title_fill = title_fill_alpha)
  dev.off()
} else {
  tiff(tiff_file, width = page_w_in, height = page_h_in, units = "in", res = 600, compression = "lzw")
  draw_box1(title_fill = title_fill_alpha)
  dev.off()
}

message("Box 1 exported: ", pdf_file, " and ", tiff_file)

##### 
##### 
##### Box 1. Rectangular lines version, italics and subscript
##### 
##### 
###############################################################################
# Box 1 | Geometric concepts and definitions used in SigPolytope
# Publication-ready export as PDF (vector) and TIFF (600 dpi)
#
# ONLY CHANGE IN THIS REVISION:
#   - Replace rounded box borders with rectangular borders
#   - NOTHING ELSE CHANGED
###############################################################################

suppressPackageStartupMessages({
  library(grid)
})

have_ragg  <- requireNamespace("ragg", quietly = TRUE)
have_cairo <- capabilities("cairo")

# ---- 0) Text sanitization -----------------------------------------------------
sanitize_ascii <- function(x) {
  x <- gsub("\u2013|\u2014", "-", x, perl = TRUE)
  x <- gsub("\u2212", "-", x, perl = TRUE)
  x <- gsub("\u00A0", " ", x, perl = TRUE)
  iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
}

# ---- 1) Content ---------------------------------------------------------------
box_title <- sanitize_ascii("Box 1 | Geometric concepts and definitions used in SigPolytope")


box_entries <- list(
  list(
    term = expression(bold("Latent space")),
    desc = "A shared multidimensional coordinate system derived from the integration of omic, phenotypic, immune, and clinical features. Each latent axis represents a composite dimension capturing correlated biological variation across these layers."
  ),
  
  # ---- NEW (semantic disambiguation anchor; no styling changes) ---------------
  list(
    term = expression(bold("Principal axes")),
    desc = "The orthogonal dimensions of the latent embedding that capture dominant sources of variance across integrated omic, phenotypic, immune, and clinical features."
  ),
  list(
    term = expression(bold("Anisotropy")),
    desc = "A geometric property of the convex hull describing unequal extension of a signature or circuitry across latent axes, reflecting directional dominance rather than positional variance."
  ),
  # ---------------------------------------------------------------------------
  
  list(
    term = expression(bold("Signature ("*italic(sig)*") and interaction ("*italic(int)*") components")),
    desc = "For each regulatory circuitry, the sig component represents the upstream regulatory entity (e.g., an ncRNA-based signature), whereas the int component represents the downstream metabolic or molecular interaction target. Each component is represented as a local polytope embedded in latent space."
  ),
  list(
    term = expression(bold("Barycenter")),
    desc = "The geometric center (mean coordinate) of a polytope, summarizing the net directional contribution of all latent dimensions for a given circuitry component."
  ),
  list(
    term = expression(bold("Barycenter distance ("*italic(d)[bary]*")")),
    desc = "The Euclidean distance between the sig and int barycenters in the shared latent space. This metric quantifies the degree of geometric separation between the regulatory and metabolic components of a circuitry."
  ),
  list(
    term = expression(bold("Barycenter-distance regimes")),
    desc = "Circuitries are discretized into four regimes based on d_bary - high concordance, moderate discordance, strong discordance, and extreme discordance. These regimes provide an interpretable categorization of continuous geometric separation and are used consistently throughout the main text and Supplementary analyses."
  ),
  list(
    term = expression(bold("Convex hull")),
    desc = "The minimal geometric envelope enclosing all vertices of a polytope. Hull geometry provides a visual and quantitative summary of the internal organization of a signature or interaction component."
  ),
  list(
    term = expression(bold("Hull volume")),
    desc = "A quantitative measure of the breadth of latent engagement, reflecting the number and diversity of latent dimensions contributing substantially to a given circuitry component."
  ),
  list(
    term = expression(bold("Symmetry ("*italic(sig)*"-"*italic(int)*" volume ratio)")),
    desc = "The relative balance between sig and int hull volumes, indicating whether latent complexity is evenly distributed between components or dominated by one side."
  ),
  list(
    term = expression(bold("Separation-with-convergence")),
    desc = "A geometric configuration in which sig and int components are spatially separated in latent space (large barycenter distance) yet exhibit concordant sign structure along specific annotated axes (e.g., survival, microenvironmental, immune). In such cases, barycenter distance captures latent-direction separation, whereas concordance reflects axis-specific functional alignment."
  )
)


for (i in seq_along(box_entries)) box_entries[[i]]$desc <- sanitize_ascii(box_entries[[i]]$desc)

# ---- 2) Layout controls -------------------------------------------------------
page_w_in <- 7
page_h_in <- 9
font_family <- "Helvetica"

border_col <- "#1F4E79"
title_fill_alpha  <- adjustcolor("#2F6FAE", alpha.f = 0.35)
title_fill_opaque <- "#2F6FAE"
body_fill          <- adjustcolor("#DCEBFA", alpha.f = 0.55)

lwd_box  <- 2.4

# centered composite box geometry (npc)
w       <- 0.88
h_body  <- 0.76
h_title <- 0.055
total_h <- h_body + h_title

x0 <- (1 - w) / 2
y_body0  <- (1 - total_h) / 2
y_title0 <- y_body0 + h_body

# padding inside body
pad_x <- 0.04
pad_y <- 0.028

# base typography (auto-fit)
fs_title <- 13.5
fs_term_base <- 10.8
fs_desc_base <- 10.5

term_lh_mult <- 1.12
desc_lh_mult <- 1.10
gap_mult     <- 0.55

target_bottom_gap_pt <- 5

# ---- 3) Wrapping --------------------------------------------------------------
wrap_to_unit_width <- function(text, max_width, gp) {
  text <- gsub("\\s+", " ", trimws(text))
  if (!nzchar(text)) return(character(0))
  
  words <- strsplit(text, " ", fixed = TRUE)[[1]]
  lines <- character(0)
  cur <- ""
  
  max_w_in <- convertWidth(max_width, "in", valueOnly = TRUE)
  for (wd in words) {
    cand <- if (cur == "") wd else paste(cur, wd)
    cand_w <- convertWidth(grobWidth(textGrob(cand, gp = gp)), "in", valueOnly = TRUE)
    if (cand_w <= max_w_in) {
      cur <- cand
    } else {
      if (nzchar(cur)) lines <- c(lines, cur)
      cur <- wd
    }
  }
  if (nzchar(cur)) lines <- c(lines, cur)
  lines
}

# ---- 4) Justification helper --------------------------------------------------
justify_line_to_width <- function(line, target_width, gp) {
  line <- gsub("\\s+", " ", trimws(line))
  if (!nzchar(line)) return(line)
  
  words <- strsplit(line, " ", fixed = TRUE)[[1]]
  if (length(words) < 2) return(line)
  
  cur_w_in <- convertWidth(grobWidth(textGrob(line, gp = gp)), "in", valueOnly = TRUE)
  tgt_w_in <- convertWidth(target_width, "in", valueOnly = TRUE)
  if (cur_w_in >= tgt_w_in) return(line)
  
  sp_w_in <- convertWidth(grobWidth(textGrob(" ", gp = gp)), "in", valueOnly = TRUE)
  if (!is.finite(sp_w_in) || sp_w_in <= 0) return(line)
  
  extra_in <- tgt_w_in - cur_w_in
  extra_spaces <- floor(extra_in / sp_w_in)
  if (extra_spaces <= 0) return(line)
  
  gaps <- length(words) - 1
  add_each <- extra_spaces %/% gaps
  rem      <- extra_spaces %% gaps
  
  out <- words[1]
  for (i in seq_len(gaps)) {
    add_i <- add_each + if (i <= rem) 1 else 0
    out <- paste0(out, paste(rep(" ", 1 + add_i), collapse = ""), words[i + 1])
  }
  out
}

# ---- 5) Token typesetting for descriptions -----------------------------------
needs_token_typeset <- function(line) {
  grepl("\\b(sig|int)\\b|\\b(d_bary|dbary)\\b", line, perl = TRUE)
}

draw_line_with_tokens <- function(line, x, y, gp_plain, gp_italic, gp_math) {
  
  m <- gregexpr("\\s+|\\S+", line, perl = TRUE)[[1]]
  tokens <- regmatches(line, list(m))[[1]]
  
  cur_x <- x
  
  for (tok in tokens) {
    
    if (grepl("^\\s+$", tok)) {
      cur_x <- cur_x + grobWidth(textGrob(tok, gp = gp_plain))
      next
    }
    
    lead  <- sub("^([\\(\\[\\{\\\"\\'\\<\\>\\.,;:!?-]*).*$", "\\1", tok, perl = TRUE)
    trail <- sub("^.*?([\\)\\]\\}\\\"\\'\\<\\>\\.,;:!?-]*)$", "\\1", tok, perl = TRUE)
    core  <- tok
    core  <- sub("^([\\(\\[\\{\\\"\\'\\<\\>\\.,;:!?-]*)", "", core, perl = TRUE)
    core  <- sub("([\\)\\]\\}\\\"\\'\\<\\>\\.,;:!?-]*)$", "", core, perl = TRUE)
    if (!nzchar(core)) core <- tok
    
    if (nzchar(lead)) {
      grid.draw(textGrob(lead, x = cur_x, y = y, just = c("left", "top"), gp = gp_plain))
      cur_x <- cur_x + grobWidth(textGrob(lead, gp = gp_plain))
    }
    
    if (core %in% c("sig", "int")) {
      grid.draw(textGrob(core, x = cur_x, y = y, just = c("left", "top"), gp = gp_italic))
      cur_x <- cur_x + grobWidth(textGrob(core, gp = gp_italic))
      
    } else if (core %in% c("d_bary", "dbary")) {
      grid.draw(textGrob(expression(italic(d)[bary]),
                         x = cur_x, y = y, just = c("left", "top"), gp = gp_math))
      cur_x <- cur_x + grobWidth(textGrob(expression(italic(d)[bary]), gp = gp_math))
      
    } else {
      grid.draw(textGrob(core, x = cur_x, y = y, just = c("left", "top"), gp = gp_plain))
      cur_x <- cur_x + grobWidth(textGrob(core, gp = gp_plain))
    }
    
    if (nzchar(trail)) {
      grid.draw(textGrob(trail, x = cur_x, y = y, just = c("left", "top"), gp = gp_plain))
      cur_x <- cur_x + grobWidth(textGrob(trail, gp = gp_plain))
    }
  }
}

# ---- 6) Build layout ----------------------------------------------------------
build_body_layout <- function(avail_w, gp_desc, lh_term_pt, lh_desc_pt, gap_pt) {
  layout <- list()
  n_terms <- 0; n_desc <- 0; n_gaps <- 0
  
  for (i in seq_along(box_entries)) {
    
    e <- box_entries[[i]]
    layout[[length(layout) + 1]] <- list(type = "term", label = e$term)
    n_terms <- n_terms + 1
    
    desc_lines <- wrap_to_unit_width(e$desc, avail_w, gp_desc)
    if (length(desc_lines) >= 1) {
      for (j in seq_along(desc_lines)) {
        is_last <- (j == length(desc_lines))
        line <- desc_lines[j]
        
        layout[[length(layout) + 1]] <- list(
          type = "desc",
          text = line,
          justify = (!is_last) && !needs_token_typeset(line),
          token_typeset = needs_token_typeset(line)
        )
        n_desc <- n_desc + 1
      }
    }
    
    if (i < length(box_entries)) {
      layout[[length(layout) + 1]] <- list(type = "gap")
      n_gaps <- n_gaps + 1
    }
  }
  
  list(
    layout = layout,
    total_height_pt = (n_terms * lh_term_pt) + (n_desc * lh_desc_pt) + (n_gaps * gap_pt)
  )
}

# ---- 7) Renderer --------------------------------------------------------------
draw_box1 <- function(title_fill) {
  
  grid.newpage()
  
  # ---- ONLY CHANGE: RECTANGULAR BOXES ----------------------------------------
  grid.rect(
    x = x0, y = y_title0, width = w, height = h_title,
    just = c("left", "bottom"),
    gp = gpar(col = border_col, fill = title_fill, lwd = lwd_box)
  )
  
  grid.rect(
    x = x0, y = y_body0, width = w, height = h_body,
    just = c("left", "bottom"),
    gp = gpar(col = border_col, fill = body_fill, lwd = lwd_box)
  )
  # ---------------------------------------------------------------------------
  
  # Title text centered
  pushViewport(viewport(x = x0, y = y_title0, width = w, height = h_title,
                        just = c("left", "bottom"), clip = "on"))
  grid.text(
    box_title,
    x = unit(0.5, "npc"), y = unit(0.5, "npc"),
    just = "center",
    gp = gpar(fontsize = fs_title, fontface = "bold", col = "black", fontfamily = font_family)
  )
  popViewport()
  
  # Body viewport
  pushViewport(viewport(
    x = x0 + pad_x, y = y_body0 + pad_y,
    width  = w - 2 * pad_x,
    height = h_body - 2 * pad_y,
    just = c("left", "bottom"),
    clip = "on"
  ))
  
  avail_w <- unit(1, "npc")
  avail_h_pt <- convertHeight(unit(1, "npc"), "pt", valueOnly = TRUE)
  
  # Auto fit-to-box scaling
  scale <- 1.0
  scale_min <- 0.80
  
  repeat {
    fs_term <- fs_term_base * scale
    fs_desc <- fs_desc_base * scale
    lh_term_pt <- fs_term * term_lh_mult
    lh_desc_pt <- fs_desc * desc_lh_mult
    gap_pt     <- fs_desc * gap_mult
    
    gp_desc_plain <- gpar(fontface = "plain", fontsize = fs_desc, col = "black", fontfamily = font_family)
    layout_obj <- build_body_layout(avail_w, gp_desc_plain, lh_term_pt, lh_desc_pt, gap_pt)
    
    if (layout_obj$total_height_pt <= (avail_h_pt - target_bottom_gap_pt)) break
    if (scale <= scale_min) break
    scale <- max(scale_min, scale * 0.985)
  }
  
  # Final GPs
  fs_term <- fs_term_base * scale
  fs_desc <- fs_desc_base * scale
  lh_term_pt <- fs_term * term_lh_mult
  lh_desc_pt <- fs_desc * desc_lh_mult
  gap_pt     <- fs_desc * gap_mult
  
  gp_term <- gpar(fontsize = fs_term, col = "black", fontfamily = font_family)
  gp_desc_plain  <- gpar(fontface = "plain",  fontsize = fs_desc, col = "black", fontfamily = font_family)
  gp_desc_italic <- gpar(fontface = "italic", fontsize = fs_desc, col = "black", fontfamily = font_family)
  gp_math <- gpar(fontface = "plain", fontsize = fs_desc, col = "black", fontfamily = font_family)
  
  layout_obj <- build_body_layout(avail_w, gp_desc_plain, lh_term_pt, lh_desc_pt, gap_pt)
  
  top_offset_pt <- max(0, avail_h_pt - layout_obj$total_height_pt - target_bottom_gap_pt)
  
  cur_x <- unit(0, "npc")
  cur_y <- unit(1, "npc") - unit(top_offset_pt, "pt")
  
  for (item in layout_obj$layout) {
    
    if (identical(item$type, "term")) {
      grid.draw(textGrob(item$label, x = cur_x, y = cur_y, just = c("left", "top"), gp = gp_term))
      cur_y <- cur_y - unit(lh_term_pt, "pt")
      
    } else if (identical(item$type, "desc")) {
      
      if (isTRUE(item$token_typeset)) {
        draw_line_with_tokens(item$text, x = cur_x, y = cur_y,
                              gp_plain = gp_desc_plain, gp_italic = gp_desc_italic, gp_math = gp_math)
      } else {
        ln <- item$text
        if (isTRUE(item$justify)) ln <- justify_line_to_width(ln, target_width = avail_w, gp = gp_desc_plain)
        grid.text(ln, x = cur_x, y = cur_y, just = c("left", "top"), gp = gp_desc_plain)
      }
      
      cur_y <- cur_y - unit(lh_desc_pt, "pt")
      
    } else if (identical(item$type, "gap")) {
      cur_y <- cur_y - unit(gap_pt, "pt")
    }
  }
  
  popViewport()
}

# ---- 8) Export ----------------------------------------------------------------
pdf_file  <- "Box1_SigPolytope_rectangular_lines.pdf"
tiff_file <- "Box1_SigPolytope_600dpi_rectangular_lines.tiff"

if (have_cairo) {
  grDevices::cairo_pdf(filename = pdf_file, width = page_w_in, height = page_h_in, family = font_family)
  draw_box1(title_fill = title_fill_alpha)
  dev.off()
} else {
  grDevices::pdf(file = pdf_file, width = page_w_in, height = page_h_in,
                 family = font_family, useDingbats = FALSE)
  draw_box1(title_fill = title_fill_opaque)
  dev.off()
}

if (have_ragg) {
  ragg::agg_tiff(filename = tiff_file, width = page_w_in, height = page_h_in,
                 units = "in", res = 600, compression = "lzw")
  draw_box1(title_fill = title_fill_alpha)
  dev.off()
} else {
  tiff(tiff_file, width = page_w_in, height = page_h_in, units = "in", res = 600, compression = "lzw")
  draw_box1(title_fill = title_fill_alpha)
  dev.off()
}

message("Box 1 exported: ", pdf_file, " and ", tiff_file)

################################################################################
# START
################################################################################



################################################################################
# Figure 2 – (V3) Conceptual Geometry (Panels A–F)
# ---------------------------------------------------------------------------
# LOCKED, REPRODUCIBLE EXPORT PIPELINE (Plotly/Kaleido → PNG → 600 dpi TIFF → PDF)
#
# REQUIREMENTS:
#   - R: plotly, geometry, magick
#   - Reticulate must be bound to Python 3.11 env with:
#       python:plotly  AND  python-kaleido==0.2.1
#   - You have already validated py_config() points to r-reticulate311 (Python 3.11)
#
# KEY CHANGE vs prior version:
#   - Removed HTML saveWidget + webshot2 browser screenshots.
#   - Uses plotly::save_image() (Kaleido) for deterministic static export.
################################################################################

suppressPackageStartupMessages({
  library(plotly)
  library(geometry)
  library(magick)
})

#-----------------------------
# 0) Optional: hard sanity check for Kaleido backend (recommended)
#-----------------------------
# If you want to enforce correctness early, uncomment:
# suppressPackageStartupMessages(library(reticulate))
# cfg <- reticulate::py_config()
# message("Python: ", cfg$python, " | Version: ", cfg$version)

#-----------------------------
# 1) Build synthetic panels A–F
#-----------------------------
set.seed(42)

generate_point_cloud <- function(type = c("compact","elongated","single_axis","multi_axis","irregular"),
                                 n = 30, scale = 1) {
  type <- match.arg(type)
  
  if (type == "compact") {
    X <- matrix(rnorm(n * 3, sd = 0.35), ncol = 3)
  } else if (type == "elongated") {
    X <- cbind(rnorm(n, sd = 1.20), rnorm(n, sd = 0.25), rnorm(n, sd = 0.25))
  } else if (type == "single_axis") {
    X <- cbind(rnorm(n, sd = 1.50), rnorm(n, sd = 0.12), rnorm(n, sd = 0.12))
  } else if (type == "multi_axis") {
    X <- cbind(rnorm(n, sd = 1.00), rnorm(n, sd = 0.95), rnorm(n, sd = 0.15))
  } else if (type == "irregular") {
    n1 <- floor(n / 3); n2 <- floor(n / 3); n3 <- n - n1 - n2
    X1 <- matrix(rnorm(n1 * 3, mean = -1.2, sd = 0.25), ncol = 3)
    X2 <- matrix(rnorm(n2 * 3, mean =  0.0, sd = 0.30), ncol = 3)
    X3 <- matrix(rnorm(n3 * 3, mean =  1.2, sd = 0.25), ncol = 3)
    X  <- rbind(X1, X2, X3)
  }
  
  X * scale
}

compute_hull <- function(X) {
  ch <- geometry::convhulln(X, options = "FA")
  list(points = X, facets = ch$hull, area = ch$area, vol = ch$vol)
}

plot_hull_3d <- function(hull_obj,
                         title = NULL,
                         point_color = "gray40",
                         point_size  = 10,
                         hull_opacity = 0.35) {
  X <- hull_obj$points
  F <- hull_obj$facets
  stopifnot(is.matrix(X), ncol(X) == 3)
  stopifnot(is.matrix(F), ncol(F) == 3)
  
  plot_ly() %>%
    add_trace(
      x = X[,1], y = X[,2], z = X[,3],
      type = "scatter3d",
      mode = "markers",
      marker = list(size = point_size, color = point_color),
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "mesh3d",
      x = X[,1], y = X[,2], z = X[,3],
      i = F[,1] - 1, j = F[,2] - 1, k = F[,3] - 1,
      opacity = hull_opacity,
      showlegend = FALSE
    ) %>%
    layout(
      title = list(text = title),
      margin = list(l = 0, r = 0, b = 0, t = 20),
      scene = list(
        xaxis = list(title = "Latent axis 1"),
        yaxis = list(title = "Latent axis 2"),
        zaxis = list(title = "Latent axis 3")
      )
    ) %>%
    config(displayModeBar = FALSE)
}

# Build data/hulls
X_A <- generate_point_cloud("compact",     n = 30)
X_B <- generate_point_cloud("compact",     n = 30)
X_C <- generate_point_cloud("elongated",   n = 30)
X_D <- generate_point_cloud("single_axis", n = 30)
X_E <- generate_point_cloud("multi_axis",  n = 30)
X_F <- generate_point_cloud("irregular",   n = 30)

h_A <- compute_hull(X_A); h_B <- compute_hull(X_B); h_C <- compute_hull(X_C)
h_D <- compute_hull(X_D); h_E <- compute_hull(X_E); h_F <- compute_hull(X_F)

fig_A <- plot_hull_3d(h_A, title = "(A) Conceptual projection of a multidimensional omic signature")
fig_B <- plot_hull_3d(h_B, title = "(B) Compact and isotropic geometric organization")
fig_C <- plot_hull_3d(h_C, title = "(C) Elongated and anisotropic geometric organization")
fig_D <- plot_hull_3d(h_D, title = "(D) Coherent single-axis regulatory organization")
fig_E <- plot_hull_3d(h_E, title = "(E) Directional multi-axis regulatory organization")
fig_F <- plot_hull_3d(h_F, title = "(F) Heterogeneous and modular organizational structure")

fig_list <- list(A = fig_A, B = fig_B, C = fig_C, D = fig_D, E = fig_E, F = fig_F)

#-----------------------------
# 2) Per-panel hull colors (LOCKED)
#-----------------------------
hull_palette <- c(
  A = "rgba(243, 156,  18, 0.35)",
  B = "rgba(231,  76,  60, 0.35)",
  C = "rgba(160, 140, 120, 0.35)",
  D = "rgba(140, 140, 140, 0.35)",
  E = "rgba( 52, 152, 219, 0.35)",
  F = "rgba(243, 156,  18, 0.35)"
)

point_color_locked <- "rgba(31, 119, 180, 1.0)"

#-----------------------------
# 3) LOCK STYLE + paper axis titles; disable 3D axis titles AND tick labels
#    + recreate conceptual tick numbers as PAPER annotations (export-stable)
#-----------------------------
lock_panel_style <- function(fig, hull_rgba,
                             point_color  = point_color_locked,
                             point_size   = 10,
                             
                             font_family  = "Arial",
                             base_font    = 26,
                             axis_title   = 50,
                             
                             tick_font    = 50,
                             title_size   = 75,
                             title_y      = 0.85,
                             top_margin_t = 20,
                             
                             xlab_y = 0.02,
                             ylab_x = 0.02,
                             zlab_x = 0.98,
                             
                             margin_l = 40,
                             margin_r = 40,
                             margin_b = 40,
                             
                             scene_domain_x = c(0.06, 0.94),
                             scene_domain_y = c(0.16, 0.92),
                             
                             hull_opacity = 0.35,
                             
                             tick_vals = c(-2, 0, 2),
                             
                             xtick_x0 = 0.26,
                             xtick_x1 = 0.74,
                             xtick_y  = 0.105,
                             
                             ytick_x  = 0.095,
                             ytick_y0 = 0.28,
                             ytick_y1 = 0.72,
                             
                             ztick_x  = 0.905,
                             ztick_y0 = 0.28,
                             ztick_y1 = 0.72,
                             
                             xtick_yshift = 0,
                             ytick_xshift = 0,
                             ztick_xshift = 0) {
  
  stopifnot(inherits(fig, "plotly"))
  
  pb <- plotly_build(fig)
  if (is.null(pb$x$data) || length(pb$x$data) == 0) {
    stop("plotly_build() produced no traces; invalid plotly object.")
  }
  
  # Extract title text, remove default title
  this_title <- ""
  if (!is.null(pb$x$layout$title) && !is.null(pb$x$layout$title$text)) {
    this_title <- pb$x$layout$title$text
  }
  pb$x$layout$title <- NULL
  
  # Edit traces
  for (i in seq_along(pb$x$data)) {
    tr <- pb$x$data[[i]]
    
    if (identical(tr$type, "scatter3d")) {
      if (is.null(tr$marker)) tr$marker <- list()
      tr$marker$color <- point_color
      tr$marker$size  <- point_size
      pb$x$data[[i]] <- tr
    }
    
    if (identical(tr$type, "mesh3d")) {
      nfaces <- length(tr$i)
      tr$facecolor <- rep(hull_rgba, nfaces)
      tr$opacity   <- hull_opacity
      
      tr$flatshading   <- TRUE
      tr$lighting      <- list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0)
      tr$lightposition <- list(x = 0, y = 0, z = 1e6)
      
      tr$mode  <- NULL
      tr$color <- NULL
      pb$x$data[[i]] <- tr
    }
  }
  
  # Panel title annotation
  ann_title <- list(
    text      = this_title,
    x         = 0.50,
    y         = title_y,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "bottom",
    showarrow = FALSE,
    font      = list(family = font_family, size = title_size, color = "black")
  )
  
  # Axis-title annotations (paper coords)
  ann_x <- list(
    text      = "Latent axis 1",
    x         = 0.50,
    y         = xlab_y,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "top",
    showarrow = FALSE,
    font      = list(family = font_family, size = axis_title, color = "black")
  )
  
  ann_y <- list(
    text      = "Latent axis 2",
    x         = ylab_x,
    y         = 0.50,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "middle",
    textangle = -90,
    showarrow = FALSE,
    font      = list(family = font_family, size = axis_title, color = "black")
  )
  
  ann_z <- list(
    text      = "Latent axis 3",
    x         = zlab_x,
    y         = 0.50,
    xref      = "paper",
    yref      = "paper",
    xanchor   = "center",
    yanchor   = "middle",
    textangle = 90,
    showarrow = FALSE,
    font      = list(family = font_family, size = axis_title, color = "black")
  )
  
  fmt_tick <- function(v) {
    if (abs(v) < 1e-12) return("0")
    if (abs(v - round(v)) < 1e-12) return(as.character(as.integer(round(v))))
    sprintf("%.1f", v)
  }
  
  x_positions <- seq(xtick_x0, xtick_x1, length.out = length(tick_vals))
  y_positions <- seq(ytick_y0, ytick_y1, length.out = length(tick_vals))
  z_positions <- seq(ztick_y0, ztick_y1, length.out = length(tick_vals))
  
  ann_xticks <- lapply(seq_along(tick_vals), function(k) {
    list(
      text      = fmt_tick(tick_vals[k]),
      x         = x_positions[k],
      y         = xtick_y,
      xref      = "paper",
      yref      = "paper",
      xanchor   = "center",
      yanchor   = "top",
      yshift    = xtick_yshift,
      showarrow = FALSE,
      font      = list(family = font_family, size = tick_font, color = "black")
    )
  })
  
  ann_yticks <- lapply(seq_along(tick_vals), function(k) {
    list(
      text      = fmt_tick(tick_vals[k]),
      x         = ytick_x,
      y         = y_positions[k],
      xref      = "paper",
      yref      = "paper",
      xanchor   = "right",
      yanchor   = "middle",
      xshift    = ytick_xshift,
      showarrow = FALSE,
      font      = list(family = font_family, size = tick_font, color = "black")
    )
  })
  
  ann_zticks <- lapply(seq_along(tick_vals), function(k) {
    list(
      text      = fmt_tick(tick_vals[k]),
      x         = ztick_x,
      y         = z_positions[k],
      xref      = "paper",
      yref      = "paper",
      xanchor   = "left",
      yanchor   = "middle",
      xshift    = ztick_xshift,
      showarrow = FALSE,
      font      = list(family = font_family, size = tick_font, color = "black")
    )
  })
  
  pb$x$layout$annotations <- c(
    list(ann_title, ann_x, ann_y, ann_z),
    ann_xticks,
    ann_yticks,
    ann_zticks
  )
  
  axis3d_no_ticks <- list(
    title = list(text = ""),
    showticklabels = FALSE,
    ticks = "",
    showticksuffix = "none",
    showexponent = "none",
    
    showgrid       = TRUE,
    gridcolor      = "rgba(0,0,0,0.35)",
    gridwidth      = 2,
    
    showbackground  = TRUE,
    backgroundcolor = "rgba(0,0,0,0.04)",
    
    showline   = TRUE,
    linecolor  = "rgba(0,0,0,0.35)",
    linewidth  = 2,
    
    zeroline      = TRUE,
    zerolinecolor = "rgba(0,0,0,0.45)",
    zerolinewidth = 2
  )
  
  pb <- pb %>%
    layout(
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      font = list(family = font_family, size = base_font, color = "black"),
      
      margin = list(l = margin_l, r = margin_r, b = margin_b, t = top_margin_t),
      
      scene = list(
        bgcolor = "white",
        domain  = list(x = scene_domain_x, y = scene_domain_y),
        xaxis   = axis3d_no_ticks,
        yaxis   = axis3d_no_ticks,
        zaxis   = axis3d_no_ticks
      )
    ) %>%
    config(displayModeBar = FALSE)
  
  pb
}

# Apply locked style A–F
fig_list_locked <- fig_list
for (nm in names(fig_list_locked)) {
  fig_list_locked[[nm]] <- lock_panel_style(
    fig       = fig_list_locked[[nm]],
    hull_rgba = hull_palette[[nm]]
  )
}

#-----------------------------
# 4) Output folders
#-----------------------------
out_png  <- "Figure2_Panels_PNG"
out_tif  <- "Figure2_Panels_TIFF_600dpi"
out_comp <- "Figure2_Composite"
dir.create(out_png,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_tif,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_comp, showWarnings = FALSE, recursive = TRUE)

#-----------------------------
# 5) Export parameters
#-----------------------------
dpi_target <- 600
panel_width_in  <- 3.5
panel_height_in <- 3.0
panel_w_px <- as.integer(round(panel_width_in  * dpi_target))
panel_h_px <- as.integer(round(panel_height_in * dpi_target))

#-----------------------------
# 6) Exporter: Plotly (Kaleido) -> PNG -> TIFF
#-----------------------------
export_plotly_panel <- function(fig, nm,
                                out_png, out_tif,
                                w_px, h_px,
                                dpi = 600) {
  
  png_file  <- file.path(out_png, sprintf("Figure2_%s.png",  nm))
  tif_file  <- file.path(out_tif, sprintf("Figure2_%s_%ddpi.tif", nm, dpi))
  
  # Deterministic static render via Kaleido (no browser/webshot2)
  plotly::save_image(fig, file = png_file, width = w_px, height = h_px, scale = 1)
  if (!file.exists(png_file)) stop("PNG not created: ", png_file)
  
  img <- magick::image_read(png_file)
  magick::image_write(
    img,
    path        = tif_file,
    format      = "tiff",
    compression = "lzw",
    density     = paste0(dpi, "x", dpi)
  )
  if (!file.exists(tif_file)) stop("TIFF not created: ", tif_file)
  
  list(png = png_file, tif = tif_file)
}

#-----------------------------
# 7) Export A–F
#-----------------------------
exports <- setNames(vector("list", length(fig_list_locked)), names(fig_list_locked))

for (nm in names(fig_list_locked)) {
  message("\n=== Exporting panel ", nm, " ===")
  exports[[nm]] <- export_plotly_panel(
    fig     = fig_list_locked[[nm]],
    nm      = nm,
    out_png = out_png,
    out_tif = out_tif,
    w_px    = panel_w_px,
    h_px    = panel_h_px,
    dpi     = dpi_target
  )
  message("OK ", nm, " -> ", exports[[nm]]$tif)
}

#-----------------------------
# 8) Composite 2×3 TIFF + PDF
#-----------------------------
panel_tifs <- vapply(exports, function(x) x$tif, FUN.VALUE = character(1))
panel_tifs <- panel_tifs[c("A","B","C","D","E","F")]

if (any(!file.exists(panel_tifs))) {
  stop("Composite aborted: not all panel TIFF files exist on disk.")
}

imgs <- lapply(panel_tifs, magick::image_read)

pad_px <- 30
imgs <- lapply(imgs, function(im) magick::image_border(im, "white", paste0(pad_px, "x", pad_px)))

row1 <- magick::image_append(magick::image_join(imgs[1:3]), stack = FALSE)
row2 <- magick::image_append(magick::image_join(imgs[4:6]), stack = FALSE)
comp <- magick::image_append(magick::image_join(list(row1, row2)), stack = TRUE)

comp_tif <- file.path(out_comp, sprintf("Figure2_Composite_2x3_%ddpi.tif", dpi_target))
magick::image_write(
  comp,
  path        = comp_tif,
  format      = "tiff",
  compression = "lzw",
  density     = paste0(dpi_target, "x", dpi_target)
)
if (!file.exists(comp_tif)) stop("Composite TIFF not created: ", comp_tif)

# PDF composite (raster embedded in PDF wrapper; acceptable for many journals)
comp_pdf_target <- file.path(out_comp, "Figure2_Composite_2x3.pdf")
magick::image_write(comp, path = comp_pdf_target, format = "pdf")
if (!file.exists(comp_pdf_target)) stop("Composite PDF not created: ", comp_pdf_target)

message("\nDONE.\nPanels written to:\n  PNG : ", normalizePath(out_png),
        "\n  TIFF: ", normalizePath(out_tif),
        "\nComposite:\n  ", normalizePath(comp_tif),
        "\n  ", normalizePath(comp_pdf_target))

################################################################################
# END
################################################################################


################################################################################
# START
################################################################################

###############################################################################
# Figure 2 – (V4) Conceptual Geometry (Panels A–F)
# NATIVE AXES + PER-PANEL COLORS + PUBLICATION-SCALE FONTS
# Deterministic export: Plotly/Kaleido -> PNG -> 600 dpi TIFF (+ PDF composite)
###############################################################################

suppressPackageStartupMessages({
  library(plotly)
  library(geometry)
  library(htmlwidgets)
  library(magick)
  library(reticulate)
})

set.seed(42)

print(reticulate::py_config())

generate_point_cloud <- function(type = c("compact","elongated","single_axis","multi_axis","irregular"),
                                 n = 30, scale = 1) {
  type <- match.arg(type)
  if (type == "compact") {
    X <- matrix(rnorm(n * 3, sd = 0.35), ncol = 3)
  } else if (type == "elongated") {
    X <- cbind(rnorm(n, sd = 1.20), rnorm(n, sd = 0.25), rnorm(n, sd = 0.25))
  } else if (type == "single_axis") {
    X <- cbind(rnorm(n, sd = 1.50), rnorm(n, sd = 0.12), rnorm(n, sd = 0.12))
  } else if (type == "multi_axis") {
    X <- cbind(rnorm(n, sd = 1.00), rnorm(n, sd = 0.95), rnorm(n, sd = 0.15))
  } else {
    n1 <- floor(n / 3); n2 <- floor(n / 3); n3 <- n - n1 - n2
    X1 <- matrix(rnorm(n1 * 3, mean = -1.2, sd = 0.25), ncol = 3)
    X2 <- matrix(rnorm(n2 * 3, mean =  0.0, sd = 0.30), ncol = 3)
    X3 <- matrix(rnorm(n3 * 3, mean =  1.2, sd = 0.25), ncol = 3)
    X  <- rbind(X1, X2, X3)
  }
  X * scale
}

compute_hull <- function(X) {
  ch <- geometry::convhulln(X, options = "FA")
  list(points = X, facets = ch$hull, area = ch$area, vol = ch$vol)
}

hull_palette <- c(
  A = "rgba(243, 156,  18, 0.35)",
  B = "rgba(231,  76,  60, 0.35)",
  C = "rgba(160, 140, 120, 0.35)",
  D = "rgba(140, 140, 140, 0.35)",
  E = "rgba( 52, 152, 219, 0.35)",
  F = "rgba(243, 156,  18, 0.35)"
)
point_color_locked <- "rgba(31, 119, 180, 1.0)"

# Publication-scale fonts (PIXELS; must be large for 600 dpi export)
FONT_FAMILY <- "Arial"
TITLE_FONT  <- 70   # panel title (A–F) inside each subplot
AXIS_TITLE  <- 60   # "Latent axis 1–3"
TICK_FONT   <- 25   # tick labels

GRID_COLOR  <- "rgba(0,0,0,0.45)"
GRID_WIDTH  <- 2
PLANE_COLOR <- "rgba(0,0,0,0.05)"
AXISLINE_COLOR <- "rgba(0,0,0,0.55)"
ZEROLINE_COLOR <- "rgba(0,0,0,0.60)"

axis3d_native <- function() {
  list(
    title = list(
      text = NULL,
      font = list(family = FONT_FAMILY, size = AXIS_TITLE, color = "black"),
      standoff = 28          # <-- KEY: pushes title away from tick labels
    ),
    tickfont = list(family = FONT_FAMILY, size = TICK_FONT, color = "black"),
    showticklabels = TRUE,
    ticks = "outside",
    ticklen = 10,
    tickwidth = 3,
    tickpadding = 14,        # <-- KEY: gives tick labels breathing room
    showgrid = TRUE,
    gridcolor = GRID_COLOR,
    gridwidth = GRID_WIDTH,
    showbackground = TRUE,
    backgroundcolor = PLANE_COLOR,
    showline = TRUE,
    linecolor = AXISLINE_COLOR,
    linewidth = 2,
    zeroline = TRUE,
    zerolinecolor = ZEROLINE_COLOR,
    zerolinewidth = 2
  )
}

plot_hull_3d <- function(hull_obj, title, hull_rgba,
                         point_color = point_color_locked,
                         point_size = 10,
                         hull_opacity = 0.35) {
  
  X <- hull_obj$points
  F <- hull_obj$facets
  stopifnot(is.matrix(X), ncol(X) == 3, is.matrix(F), ncol(F) == 3)
  
  plot_ly() %>%
    add_trace(
      x = X[,1], y = X[,2], z = X[,3],
      type = "scatter3d",
      mode = "markers",
      marker = list(size = point_size, color = point_color),
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "mesh3d",
      x = X[,1], y = X[,2], z = X[,3],
      i = F[,1] - 1, j = F[,2] - 1, k = F[,3] - 1,
      facecolor = rep(hull_rgba, nrow(F)),
      opacity = hull_opacity,
      flatshading = TRUE,
      lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
      lightposition = list(x = 0, y = 0, z = 1e6),
      showlegend = FALSE
    ) %>%
    layout(
      title = list(
        text = title,
        x = 0.02, xanchor = "left",
        y = 0.975, yanchor = "top",   # <-- KEY: prevents clipping at the raster edge
        font = list(family = FONT_FAMILY, size = TITLE_FONT, color = "black")
      ),
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      margin = list(l = 0, r = 0, b = 0, t = 110),
      scene = list(
        bgcolor = "white",
        xaxis = modifyList(axis3d_native(), list(title = list(text = "Latent axis 1"))),
        yaxis = modifyList(axis3d_native(), list(title = list(text = "Latent axis 2"))),
        zaxis = modifyList(axis3d_native(), list(title = list(text = "Latent axis 3")))
      )
    ) %>%
    config(displayModeBar = FALSE)
}

X_A <- generate_point_cloud("compact",     n = 30)
X_B <- generate_point_cloud("compact",     n = 30)
X_C <- generate_point_cloud("elongated",   n = 30)
X_D <- generate_point_cloud("single_axis", n = 30)
X_E <- generate_point_cloud("multi_axis",  n = 30)
X_F <- generate_point_cloud("irregular",   n = 30)

h_A <- compute_hull(X_A); h_B <- compute_hull(X_B); h_C <- compute_hull(X_C)
h_D <- compute_hull(X_D); h_E <- compute_hull(X_E); h_F <- compute_hull(X_F)

fig_list <- list(
  A = plot_hull_3d(h_A, "(A) Conceptual projection of a multidimensional omic signature", hull_palette["A"]),
  B = plot_hull_3d(h_B, "(B) Compact and isotropic geometric organization",             hull_palette["B"]),
  C = plot_hull_3d(h_C, "(C) Elongated and anisotropic geometric organization",         hull_palette["C"]),
  D = plot_hull_3d(h_D, "(D) Coherent single-axis regulatory organization",             hull_palette["D"]),
  E = plot_hull_3d(h_E, "(E) Directional multi-axis regulatory organization",           hull_palette["E"]),
  F = plot_hull_3d(h_F, "(F) Heterogeneous and modular organizational structure",       hull_palette["F"])
)

dpi_target <- 600
panel_width_in  <- 3.5
panel_height_in <- 3.0
panel_w_px <- as.integer(round(panel_width_in  * dpi_target))
panel_h_px <- as.integer(round(panel_height_in * dpi_target))
kaleido_scale <- 2

out_html <- "Figure2_nativeAxes_HTML"
out_png  <- "Figure2_nativeAxes_PNG"
out_tif  <- "Figure2_nativeAxes_TIFF_600dpi"
out_comp <- "Figure2_nativeAxes_Composite"
dir.create(out_html, showWarnings = FALSE, recursive = TRUE)
dir.create(out_png,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_tif,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_comp, showWarnings = FALSE, recursive = TRUE)

export_panel <- function(fig, nm) {
  html_file <- file.path(out_html, sprintf("Figure2_%s.html", nm))
  png_file  <- file.path(out_png,  sprintf("Figure2_%s_%ddpi.png", nm, dpi_target))
  tif_file  <- file.path(out_tif,  sprintf("Figure2_%s_%ddpi.tif", nm, dpi_target))
  
  htmlwidgets::saveWidget(fig, file = html_file, selfcontained = TRUE)
  
  # Ensure fully built object (robustness)
  fig_built <- plotly_build(fig)
  
  # IMPORTANT: positional arguments (p, file) — NOT fig=
  plotly::save_image(
    fig_built,
    png_file,
    width  = panel_w_px,
    height = panel_h_px,
    scale  = kaleido_scale
  )
  
  img <- magick::image_read(png_file)
  magick::image_write(
    img,
    path = tif_file,
    format = "tiff",
    compression = "lzw",
    density = paste0(dpi_target, "x", dpi_target)
  )
  
  list(html = html_file, png = png_file, tif = tif_file)
}

exports <- list()
for (nm in names(fig_list)) {
  message("Exporting panel ", nm, " ...")
  exports[[nm]] <- export_panel(fig_list[[nm]], nm)
  message("  OK: ", exports[[nm]]$tif)
}

panel_tifs <- vapply(exports[c("A","B","C","D","E","F")], `[[`, FUN.VALUE = character(1), "tif")
imgs <- lapply(panel_tifs, magick::image_read)

pad_px <- 30
imgs <- lapply(imgs, function(im) magick::image_border(im, "white", paste0(pad_px, "x", pad_px)))

row1 <- magick::image_append(magick::image_join(imgs[1:3]), stack = FALSE)
row2 <- magick::image_append(magick::image_join(imgs[4:6]), stack = FALSE)
comp <- magick::image_append(magick::image_join(list(row1, row2)), stack = TRUE)

comp_tif <- file.path(out_comp, sprintf("Figure2_Composite_2x3_%ddpi.tif", dpi_target))
magick::image_write(comp, path = comp_tif, format = "tiff", compression = "lzw",
                    density = paste0(dpi_target, "x", dpi_target))

comp_pdf <- file.path(out_comp, "Figure2_Composite_2x3_nativeAxesColored.pdf")
magick::image_write(comp, path = comp_pdf, format = "pdf")

message("\nDONE.\nPanels:\n  ", normalizePath(out_tif), "\nComposite:\n  ",
        normalizePath(comp_tif), "\n  ", normalizePath(comp_pdf))

################################################################################
# END
################################################################################

################################################################################
# START
################################################################################

###############################################################################
# Figure 2 – V6 almost final - (HUGE sized) Conceptual Geometry (Panels A–F)
# NATIVE AXES (TRUE) + PER-PANEL COLORS + LARGE EXPORT-SAFE FONTS
# Deterministic export: Plotly/Kaleido -> PNG -> 600 dpi TIFF (+ PDF composite)
###############################################################################

suppressPackageStartupMessages({
  library(plotly)
  library(geometry)
  library(htmlwidgets)
  library(magick)
  library(reticulate)
})

options(timeout = 600)  # seconds (10 minutes)

set.seed(42)
print(reticulate::py_config())

# -----------------------------
# Synthetic point clouds
# -----------------------------
generate_point_cloud <- function(type = c("compact","elongated","single_axis","multi_axis","irregular"),
                                 n = 30, scale = 1) {
  type <- match.arg(type)
  if (type == "compact") {
    X <- matrix(rnorm(n * 3, sd = 0.35), ncol = 3)
  } else if (type == "elongated") {
    X <- cbind(rnorm(n, sd = 1.20), rnorm(n, sd = 0.25), rnorm(n, sd = 0.25))
  } else if (type == "single_axis") {
    X <- cbind(rnorm(n, sd = 1.50), rnorm(n, sd = 0.12), rnorm(n, sd = 0.12))
  } else if (type == "multi_axis") {
    X <- cbind(rnorm(n, sd = 1.00), rnorm(n, sd = 0.95), rnorm(n, sd = 0.15))
  } else {
    n1 <- floor(n / 3); n2 <- floor(n / 3); n3 <- n - n1 - n2
    X1 <- matrix(rnorm(n1 * 3, mean = -1.2, sd = 0.25), ncol = 3)
    X2 <- matrix(rnorm(n2 * 3, mean =  0.0, sd = 0.30), ncol = 3)
    X3 <- matrix(rnorm(n3 * 3, mean =  1.2, sd = 0.25), ncol = 3)
    X  <- rbind(X1, X2, X3)
  }
  X * scale
}

compute_hull <- function(X) {
  ch <- geometry::convhulln(X, options = "FA")
  list(points = X, facets = ch$hull, area = ch$area, vol = ch$vol)
}

# -----------------------------
# Styling (native axes)
# -----------------------------
hull_palette <- c(
  A = "rgba(243, 156,  18, 0.35)",
  B = "rgba(231,  76,  60, 0.35)",
  C = "rgba(160, 140, 120, 0.35)",
  D = "rgba(140, 140, 140, 0.35)",
  E = "rgba( 52, 152, 219, 0.35)",
  F = "rgba(243, 156,  18, 0.35)"
)
point_color_locked <- "rgba(31, 119, 180, 1.0)"

FONT_FAMILY <- "Arial"

# These must be LARGE in Plotly/Kaleido exports (they are px, not "pt")
PANEL_TITLE_FONT <- 70
AXIS_TITLE_FONT  <- 40
TICK_FONT        <- 22   # <-- you noted 25 was still too small; keep high if needed

GRID_COLOR       <- "rgba(0,0,0,0.45)"
GRID_WIDTH       <- 2
PLANE_COLOR      <- "rgba(0,0,0,0.05)"
AXISLINE_COLOR   <- "rgba(0,0,0,0.55)"
ZEROLINE_COLOR   <- "rgba(0,0,0,0.60)"

# IMPORTANT: Build the FULL axis object with title text included.
# Do NOT overwrite title later (avoid modifyList(title=...)).
axis3d_native <- function(title_text) {
  list(
    title = list(
      text = title_text,
      font = list(family = FONT_FAMILY, size = AXIS_TITLE_FONT, color = "black"),
      standoff = 36               # separation title <-> ticks
    ),
    tickfont = list(family = FONT_FAMILY, size = TICK_FONT, color = "black"),
    showticklabels = TRUE,
    ticks = "outside",
    ticklen = 10,
    tickwidth = 3,
    tickpadding = 18,            # breathing room for tick labels
    showgrid = TRUE,
    gridcolor = GRID_COLOR,
    gridwidth = GRID_WIDTH,
    showbackground = TRUE,
    backgroundcolor = PLANE_COLOR,
    showline = TRUE,
    linecolor = AXISLINE_COLOR,
    linewidth = 2,
    zeroline = TRUE,
    zerolinecolor = ZEROLINE_COLOR,
    zerolinewidth = 2
  )
}

plot_hull_3d <- function(hull_obj, panel_title, hull_rgba,
                         point_color = point_color_locked,
                         point_size = 10,
                         hull_opacity = 0.35) {
  
  X <- hull_obj$points
  F <- hull_obj$facets
  stopifnot(is.matrix(X), ncol(X) == 3, is.matrix(F), ncol(F) == 3)
  
  plot_ly() %>%
    add_trace(
      x = X[,1], y = X[,2], z = X[,3],
      type = "scatter3d", mode = "markers",
      marker = list(size = point_size, color = point_color),
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "mesh3d",
      x = X[,1], y = X[,2], z = X[,3],
      i = F[,1] - 1, j = F[,2] - 1, k = F[,3] - 1,
      facecolor = rep(hull_rgba, nrow(F)),
      opacity = hull_opacity,
      flatshading = TRUE,
      lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
      lightposition = list(x = 0, y = 0, z = 1e6),
      showlegend = FALSE
    ) %>%
    layout(
      title = list(
        text = panel_title,
        x = 0.02, xanchor = "left",
        y = 0.98, yanchor = "top",
        font = list(family = FONT_FAMILY, size = PANEL_TITLE_FONT, color = "black")
      ),
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      margin = list(l = 0, r = 0, b = 0, t = 90),
      scene = list(
        bgcolor = "white",
        xaxis = axis3d_native("Latent axis 1"),
        yaxis = axis3d_native("Latent axis 2"),
        zaxis = axis3d_native("Latent axis 3")
      )
    ) %>%
    config(displayModeBar = FALSE)
}

# -----------------------------
# Build panels
# -----------------------------
X_A <- generate_point_cloud("compact",     n = 30)
X_B <- generate_point_cloud("compact",     n = 30)
X_C <- generate_point_cloud("elongated",   n = 30)
X_D <- generate_point_cloud("single_axis", n = 30)
X_E <- generate_point_cloud("multi_axis",  n = 30)
X_F <- generate_point_cloud("irregular",   n = 30)

h_A <- compute_hull(X_A); h_B <- compute_hull(X_B); h_C <- compute_hull(X_C)
h_D <- compute_hull(X_D); h_E <- compute_hull(X_E); h_F <- compute_hull(X_F)

fig_list <- list(
  A = plot_hull_3d(h_A, "(A) Conceptual projection of a multidimensional omic signature", hull_palette["A"]),
  B = plot_hull_3d(h_B, "(B) Compact and isotropic geometric organization",               hull_palette["B"]),
  C = plot_hull_3d(h_C, "(C) Elongated and anisotropic geometric organization",           hull_palette["C"]),
  D = plot_hull_3d(h_D, "(D) Coherent single-axis regulatory organization",               hull_palette["D"]),
  E = plot_hull_3d(h_E, "(E) Directional multi-axis regulatory organization",             hull_palette["E"]),
  F = plot_hull_3d(h_F, "(F) Heterogeneous and modular organizational structure",         hull_palette["F"])
)

# -----------------------------
# Export (Kaleido)
# -----------------------------
dpi_target <- 600
panel_width_in  <- 3.5
panel_height_in <- 3.0
panel_w_px <- as.integer(round(panel_width_in  * dpi_target))
panel_h_px <- as.integer(round(panel_height_in * dpi_target))
kaleido_scale <- 2  # increase to 3 if you want heavier raster density

out_html <- "Figure2_nativeAxes_HTML"
out_png  <- "Figure2_nativeAxes_PNG"
out_tif  <- "Figure2_nativeAxes_TIFF_600dpi"
out_comp <- "Figure2_nativeAxes_Composite"
dir.create(out_html, showWarnings = FALSE, recursive = TRUE)
dir.create(out_png,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_tif,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_comp, showWarnings = FALSE, recursive = TRUE)

export_panel <- function(fig, nm) {
  html_file <- file.path(out_html, sprintf("Figure2_%s.html", nm))
  png_file  <- file.path(out_png,  sprintf("Figure2_%s_%ddpi.png", nm, dpi_target))
  tif_file  <- file.path(out_tif,  sprintf("Figure2_%s_%ddpi.tif", nm, dpi_target))
  
  htmlwidgets::saveWidget(fig, file = html_file, selfcontained = TRUE)
  
  fig_built <- plotly_build(fig)
  plotly::save_image(
    fig_built,
    png_file,
    width  = panel_w_px,
    height = panel_h_px,
    scale  = kaleido_scale
  )
  
  img <- magick::image_read(png_file)
  magick::image_write(
    img, path = tif_file,
    format = "tiff",
    compression = "lzw",
    density = paste0(dpi_target, "x", dpi_target)
  )
  
  list(html = html_file, png = png_file, tif = tif_file)
}

exports <- list()
for (nm in names(fig_list)) {
  message("Exporting panel ", nm, " ...")
  exports[[nm]] <- export_panel(fig_list[[nm]], nm)
  message("  OK: ", exports[[nm]]$tif)
}

# -----------------------------
# Composite 2×3 WITH TOP + MID + BOTTOM GAPS (publication padding)
# -----------------------------
panel_tifs <- vapply(exports[c("A","B","C","D","E","F")], `[[`, FUN.VALUE = character(1), "tif")
imgs <- lapply(panel_tifs, magick::image_read)

pad_px <- 30
imgs <- lapply(imgs, function(im) magick::image_border(im, "white", paste0(pad_px, "x", pad_px)))

row1 <- magick::image_append(magick::image_join(imgs[1:3]), stack = FALSE)
row2 <- magick::image_append(magick::image_join(imgs[4:6]), stack = FALSE)

# Force identical widths for clean spacer insertion
w1 <- magick::image_info(row1)$width
w2 <- magick::image_info(row2)$width
wmax <- max(w1, w2)

row1 <- magick::image_extent(
  row1,
  geometry = paste0(wmax, "x", magick::image_info(row1)$height),
  gravity  = "center",
  color    = "white"
)
row2 <- magick::image_extent(
  row2,
  geometry = paste0(wmax, "x", magick::image_info(row2)$height),
  gravity  = "center",
  color    = "white"
)

# >>> USER-TUNABLE GAPS (px) <<<
TOP_GAP_PX <- 140   # adds headroom above the upper-row titles
MID_GAP_PX <- 160   # your working midline spacing
BOT_GAP_PX <- 140   # adds footroom below the lower-row axes/ticks

top_gap <- magick::image_blank(width = wmax, height = TOP_GAP_PX, color = "white")
mid_gap <- magick::image_blank(width = wmax, height = MID_GAP_PX, color = "white")
bot_gap <- magick::image_blank(width = wmax, height = BOT_GAP_PX, color = "white")

comp <- magick::image_append(
  magick::image_join(list(top_gap, row1, mid_gap, row2, bot_gap)),
  stack = TRUE
)

comp_tif <- file.path(out_comp, sprintf("Figure2_Composite_2x3_nativeAxesColored.tif", dpi_target))
magick::image_write(
  comp, path = comp_tif,
  format = "tiff",
  compression = "lzw",
  density = paste0(dpi_target, "x", dpi_target)
)

comp_pdf <- file.path(out_comp, "Figure2_Composite_2x3_nativeAxesColored.pdf")
magick::image_write(comp, path = comp_pdf, format = "pdf")

message("\nDONE.\nPanels:\n  ", normalizePath(out_tif), "\nComposite:\n  ",
        normalizePath(comp_tif), "\n  ", normalizePath(comp_pdf))

################################################################################
# END
################################################################################


###############################################################################
# Setting a dedicated directory for magick temporary files 
# Rrsetting .Renviron correctly
# 
# Minimal forensic check: magick/ImageMagick temp & cache redirection (Windows)
#
# Goal:
#   Confirm that R is using the expected startup directory (~) and that the
#   MAGICK + TEMP/TMP/TMPDIR variables are active (pointing to D:/magick_tmp),
#   preventing ImageMagick pixel-cache failures ("No space left on device").
#
# Notes:
#   - .Renviron is read only at R startup. If values are missing, restart R.
#   - .Rprofile may override environment variables after .Renviron is read.
###############################################################################

## 1) Identify the effective R home directory (~) used in this session
home <- path.expand("~")
message("R home (~): ", home)

## 2) Read the .Renviron that R would use by default (if present)
renv_path <- file.path(home, ".Renviron")
message("Active .Renviron path: ", renv_path)
message("Active .Renviron exists: ", file.exists(renv_path))
if (file.exists(renv_path)) {
  cat("\n---- Contents of active .Renviron ----\n")
  cat(readLines(renv_path, warn = FALSE), sep = "\n")
  cat("\n-------------------------------------\n")
}

## 3) Check whether explicit override pointers are in use (usually empty)
cat("\n---- R env override pointers ----\n")
print(Sys.getenv(c("R_ENVIRON", "R_ENVIRON_USER")))
cat("---------------------------------\n")

## 4) Ensure the target cache/temp directory exists
dir.create("D:/magick_tmp", recursive = TRUE, showWarnings = FALSE)

## 5) Verify the runtime environment variables that matter for magick/ImageMagick
cat("\n---- Runtime env vars (must reflect your intended configuration) ----\n")
print(Sys.getenv(c(
  "MAGICK_TEMPORARY_PATH","MAGICK_TMPDIR",
  "TMPDIR","TEMP","TMP",
  "MAGICK_MEMORY_LIMIT","MAGICK_DISK_LIMIT"
)))
cat("--------------------------------------------------------------------\n")

## 6) Optional: detect presence of .Rprofile (common override point)
rprof_path <- file.path(home, ".Rprofile")
message("Active .Rprofile exists: ", file.exists(rprof_path))
if (file.exists(rprof_path)) {
  cat("\n---- Contents of active .Rprofile ----\n")
  cat(readLines(rprof_path, warn = FALSE), sep = "\n")
  cat("\n-------------------------------------\n")
}

## 7) Optional: controlled magick cache stress test (validates disk-backed cache)
##    Run only if you want an explicit end-to-end confirmation.
# suppressPackageStartupMessages(library(magick))
# img <- image_blank(width = 8000, height = 8000, color = "white")
# outfile <- file.path(tempdir(), "magick_cache_test.tiff")
# image_write(img, path = outfile, format = "tiff")
# message("Wrote: ", outfile)
# print(file.info(outfile)$size)


###############################################################################
# Figure 2 – V7 - LOCKED FINAL (EXPORT-STABLE framing; bottom clipping prevented; auto-range padding (halo) reduced)
# Conceptual Geometry (Panels A–F) — native 3D axes + per-panel hull colors
#
# Deterministic workflow:
#   Plotly (Kaleido) -> per-panel PNG -> 600 dpi LZW TIFF -> 2×3 composite (TIFF + PDF)
#
# Key stabilization decisions (final):
#   1) No scene-domain shrink (scene uses full domain):
#        SCENE_DOMAIN_X/Y = c(0, 1)
#      This avoids global “scene shrinkage” introduced by domain padding.
#
#   2) Deterministic camera with mild zoom-in:
#        CAMERA_EYE set symmetrically (x = y = z) to preserve a constant viewpoint
#        while reducing unused internal space around the scene.
#
#   3) Halo control without altering aspectmode:
#        Per-panel axis ranges are explicitly set from the data (tight_ranges_3d),
#        using a small fractional pad plus a minimum pad to avoid degenerate ranges.
#        This reduces Plotly “nice-range” expansion while preserving the default
#        aspect scaling policy (no forced 'data' or 'cube' aspectmode).
#
#   4) Raster-edge safety:
#        A small bottom margin buffer is retained to prevent Kaleido cropping at
#        the lower edge (grid/ticks) without using scene-domain padding.
#
#   5) Scope-correct cleanup:
#        Temporary 'img' objects exist only inside export_panel(); cleanup is done
#        there, avoiding warnings from removing non-existent objects at script end.
# LOCKED PARAMETERS (do not change without regenerating all panels):
#   CAMERA_EYE = (1.39, 1.39, 1.39)
#   axis ranges: tight_ranges_3d(pad_frac = 0.03, pad_min = 1e-3)
#   margins: l=0, r=0, b=10, t=10
#   composite padding: pad_px=15; TOP/MID/BOT gaps = 10/10/10
###############################################################################

suppressPackageStartupMessages({
  library(plotly)
  library(geometry)
  library(htmlwidgets)
  library(magick)
  library(reticulate)
})

options(timeout = 600)  # seconds (10 minutes)

set.seed(42)
print(reticulate::py_config())

# -----------------------------
# Synthetic point clouds
# -----------------------------
generate_point_cloud <- function(type = c("compact","elongated","single_axis","multi_axis","irregular"),
                                 n = 30, scale = 1) {
  type <- match.arg(type)
  if (type == "compact") {
    X <- matrix(rnorm(n * 3, sd = 0.35), ncol = 3)
  } else if (type == "elongated") {
    X <- cbind(rnorm(n, sd = 1.20), rnorm(n, sd = 0.25), rnorm(n, sd = 0.25))
  } else if (type == "single_axis") {
    X <- cbind(rnorm(n, sd = 1.50), rnorm(n, sd = 0.12), rnorm(n, sd = 0.12))
  } else if (type == "multi_axis") {
    X <- cbind(rnorm(n, sd = 1.00), rnorm(n, sd = 0.95), rnorm(n, sd = 0.15))
  } else {
    n1 <- floor(n / 3); n2 <- floor(n / 3); n3 <- n - n1 - n2
    X1 <- matrix(rnorm(n1 * 3, mean = -1.2, sd = 0.25), ncol = 3)
    X2 <- matrix(rnorm(n2 * 3, mean =  0.0, sd = 0.30), ncol = 3)
    X3 <- matrix(rnorm(n3 * 3, mean =  1.2, sd = 0.25), ncol = 3)
    X  <- rbind(X1, X2, X3)
  }
  X * scale
}

compute_hull <- function(X) {
  ch <- geometry::convhulln(X, options = "FA")
  list(points = X, facets = ch$hull, area = ch$area, vol = ch$vol)
}

# -----------------------------
# Styling (native axes)
# -----------------------------
hull_palette <- c(
  A = "rgba(243, 156,  18, 0.35)",
  B = "rgba(231,  76,  60, 0.35)",
  C = "rgba(160, 140, 120, 0.35)",
  D = "rgba(140, 140, 140, 0.35)",
  E = "rgba( 52, 152, 219, 0.35)",
  F = "rgba(243, 156,  18, 0.35)"
)
point_color_locked <- "rgba(31, 119, 180, 1.0)"

FONT_FAMILY <- "Arial"

PANEL_TITLE_FONT <- 70
AXIS_TITLE_FONT  <- 40
TICK_FONT        <- 22

GRID_COLOR       <- "rgba(0,0,0,0.45)"
GRID_WIDTH       <- 2
PLANE_COLOR      <- "rgba(0,0,0,0.05)"
AXISLINE_COLOR   <- "rgba(0,0,0,0.55)"
ZEROLINE_COLOR   <- "rgba(0,0,0,0.60)"

# -----------------------------
# Scene framing controls
# -----------------------------
# Scene domain uses full occupancy (no scene-domain padding) to avoid global scene shrink.
# Clipping control is handled by a small bottom margin plus conservative axis-range padding.

SCENE_DOMAIN_X <- c(0.00, 1.00)
SCENE_DOMAIN_Y <- c(0.00, 1.00)

# Lock camera for deterministic, export-stable framing (unchanged)
CAMERA_EYE <- list(x = 1.39, y = 1.39, z = 1.39)
CAMERA_UP  <- list(x = 0, y = 0, z = 1)

axis3d_native <- function(title_text) {
  list(
    title = list(
      text = title_text,
      font = list(family = FONT_FAMILY, size = AXIS_TITLE_FONT, color = "black"),
      standoff = 22
    ),
    tickfont = list(family = FONT_FAMILY, size = TICK_FONT, color = "black"),
    showticklabels = TRUE,
    ticks = "outside",
    ticklen = 7,
    tickwidth = 2,
    tickpadding = 8,
    showgrid = TRUE,
    gridcolor = GRID_COLOR,
    gridwidth = GRID_WIDTH,
    showbackground = TRUE,
    backgroundcolor = PLANE_COLOR,
    showline = TRUE,
    linecolor = AXISLINE_COLOR,
    linewidth = 2,
    zeroline = TRUE,
    zerolinecolor = ZEROLINE_COLOR,
    zerolinewidth = 2
  )
}

tight_ranges_3d <- function(X, pad_frac = 0.03, pad_min = 1e-6) {
  stopifnot(is.matrix(X), ncol(X) == 3)
  rng  <- apply(X, 2, range)
  span <- rng[2, ] - rng[1, ]
  pad  <- pmax(span * pad_frac, pad_min)
  list(
    x = c(rng[1,1] - pad[1], rng[2,1] + pad[1]),
    y = c(rng[1,2] - pad[2], rng[2,2] + pad[2]),
    z = c(rng[1,3] - pad[3], rng[2,3] + pad[3])
  )
}


plot_hull_3d <- function(hull_obj, panel_title, hull_rgba,
                         point_color = point_color_locked,
                         point_size = 10,
                         hull_opacity = 0.35) {
  
  X <- hull_obj$points
  rng <- tight_ranges_3d(X, pad_frac = 0.03, pad_min = 1e-3)
  
  F <- hull_obj$facets
  stopifnot(is.matrix(X), ncol(X) == 3, is.matrix(F), ncol(F) == 3)
  
  plot_ly() %>%
    add_trace(
      x = X[,1], y = X[,2], z = X[,3],
      type = "scatter3d", mode = "markers",
      marker = list(size = point_size, color = point_color),
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "mesh3d",
      x = X[,1], y = X[,2], z = X[,3],
      i = F[,1] - 1, j = F[,2] - 1, k = F[,3] - 1,
      facecolor = rep(hull_rgba, nrow(F)),
      opacity = hull_opacity,
      flatshading = TRUE,
      lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
      lightposition = list(x = 0, y = 0, z = 1e6),
      showlegend = FALSE
    ) %>%
    layout(
      title = list(
        text = panel_title,
        x = 0.02, xanchor = "left",
        y = 0.98, yanchor = "top",
        font = list(family = FONT_FAMILY, size = PANEL_TITLE_FONT, color = "black")
      ),
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      
      # Small raster-edge buffer (bottom margin) to prevent Kaleido cropping of the lower
      # grid/ticks while keeping the 3D scene at full domain occupancy.
      margin = list(l = 0, r = 0, b = 10, t = 10),
      
      scene = list(
        bgcolor = "white",
        
        domain = list(x = SCENE_DOMAIN_X, y = SCENE_DOMAIN_Y),
        
        camera = list(
          eye = CAMERA_EYE,
          up  = CAMERA_UP
        ),
        
        xaxis = modifyList(axis3d_native("Latent axis 1"), list(range = rng$x)),
        yaxis = modifyList(axis3d_native("Latent axis 2"), list(range = rng$y)),
        zaxis = modifyList(axis3d_native("Latent axis 3"), list(range = rng$z))
      )
    ) %>%
    config(displayModeBar = FALSE)
}

# -----------------------------
# Build panels
# -----------------------------
X_A <- generate_point_cloud("compact",     n = 30)
X_B <- generate_point_cloud("compact",     n = 30)
X_C <- generate_point_cloud("elongated",   n = 30)
X_D <- generate_point_cloud("single_axis", n = 30)
X_E <- generate_point_cloud("multi_axis",  n = 30)
X_F <- generate_point_cloud("irregular",   n = 30)

h_A <- compute_hull(X_A); h_B <- compute_hull(X_B); h_C <- compute_hull(X_C)
h_D <- compute_hull(X_D); h_E <- compute_hull(X_E); h_F <- compute_hull(X_F)

fig_list <- list(
  A = plot_hull_3d(h_A, "(A) Conceptual projection of a multidimensional omic signature", hull_palette["A"]),
  B = plot_hull_3d(h_B, "(B) Compact and isotropic geometric organization",               hull_palette["B"]),
  C = plot_hull_3d(h_C, "(C) Elongated and anisotropic geometric organization",           hull_palette["C"]),
  D = plot_hull_3d(h_D, "(D) Coherent single-axis regulatory organization",               hull_palette["D"]),
  E = plot_hull_3d(h_E, "(E) Directional multi-axis regulatory organization",             hull_palette["E"]),
  F = plot_hull_3d(h_F, "(F) Heterogeneous and modular organizational structure",         hull_palette["F"])
)

# -----------------------------
# Export (Kaleido)
# -----------------------------
dpi_target <- 600
panel_width_in  <- 3.5
panel_height_in <- 3.0
panel_w_px <- as.integer(round(panel_width_in  * dpi_target))
panel_h_px <- as.integer(round(panel_height_in * dpi_target))
kaleido_scale <- 2

out_html <- "Figure2_nativeAxes_HTML"
out_png  <- "Figure2_nativeAxes_PNG"
out_tif  <- "Figure2_nativeAxes_TIFF_600dpi"
out_comp <- "Figure2_nativeAxes_Composite"
dir.create(out_html, showWarnings = FALSE, recursive = TRUE)
dir.create(out_png,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_tif,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_comp, showWarnings = FALSE, recursive = TRUE)

export_panel <- function(fig, nm) {
  html_file <- file.path(out_html, sprintf("Figure2_%s.html", nm))
  png_file  <- file.path(out_png,  sprintf("Figure2_%s_%ddpi.png", nm, dpi_target))
  tif_file  <- file.path(out_tif,  sprintf("Figure2_%s_%ddpi.tif", nm, dpi_target))
  
  htmlwidgets::saveWidget(fig, file = html_file, selfcontained = TRUE)
  
  fig_built <- plotly_build(fig)
  plotly::save_image(
    fig_built,
    png_file,
    width  = panel_w_px,
    height = panel_h_px,
    scale  = kaleido_scale
  )
  
  img <- magick::image_read(png_file)
  magick::image_write(
    img, path = tif_file,
    format = "tiff",
    compression = "lzw",
    density = paste0(dpi_target, "x", dpi_target)
  )
  
  invisible(gc())
  
  list(html = html_file, png = png_file, tif = tif_file)
}

exports <- list()
for (nm in names(fig_list)) {
  message("Exporting panel ", nm, " ...")
  exports[[nm]] <- export_panel(fig_list[[nm]], nm)
  message("  OK: ", exports[[nm]]$tif)
}

# -----------------------------
# Composite 2×3 WITH TOP + MID + BOTTOM GAPS (publication padding)
# -----------------------------
panel_tifs <- vapply(exports[c("A","B","C","D","E","F")], `[[`, FUN.VALUE = character(1), "tif")
imgs <- lapply(panel_tifs, magick::image_read)

pad_px <- 15
imgs <- lapply(imgs, function(im) magick::image_border(im, "white", paste0(pad_px, "x", pad_px)))

row1 <- magick::image_append(magick::image_join(imgs[1:3]), stack = FALSE)
row2 <- magick::image_append(magick::image_join(imgs[4:6]), stack = FALSE)

w1 <- magick::image_info(row1)$width
w2 <- magick::image_info(row2)$width
wmax <- max(w1, w2)

row1 <- magick::image_extent(
  row1,
  geometry = paste0(wmax, "x", magick::image_info(row1)$height),
  gravity  = "center",
  color    = "white"
)
row2 <- magick::image_extent(
  row2,
  geometry = paste0(wmax, "x", magick::image_info(row2)$height),
  gravity  = "center",
  color    = "white"
)

TOP_GAP_PX <- 10
MID_GAP_PX <- 10
BOT_GAP_PX <- 10

top_gap <- magick::image_blank(width = wmax, height = TOP_GAP_PX, color = "white")
mid_gap <- magick::image_blank(width = wmax, height = MID_GAP_PX, color = "white")
bot_gap <- magick::image_blank(width = wmax, height = BOT_GAP_PX, color = "white")

comp <- magick::image_append(
  magick::image_join(list(top_gap, row1, mid_gap, row2, bot_gap)),
  stack = TRUE
)

comp_tif <- file.path(out_comp, "Figure2_Composite_2x3_nativeAxesColored.tif")
magick::image_write(
  comp, path = comp_tif,
  format = "tiff",
  compression = "lzw",
  density = paste0(dpi_target, "x", dpi_target)
)

comp_pdf <- file.path(out_comp, "Figure2_Composite_2x3_nativeAxesColored.pdf")
magick::image_write(comp, path = comp_pdf, format = "pdf")

message("\nDONE.\nPanels:\n  ", normalizePath(out_tif), "\nComposite:\n  ",
        normalizePath(comp_tif), "\n  ", normalizePath(comp_pdf))

invisible(gc())
################################################################################
# END
################################################################################


###############################################################################
# Figure 3 HTML — FINAL! (Locked)
# Baseline plotting function (universal) with a dedicated RIGHT legend column
#
# RIGHT legend column layout (paper coordinates; single column):
#   (TOP)    Plotly trace legend (colored items: vertices / hulls / barycenters)
#   (MIDDLE) Vertex/barycenter legend (boxed annotation; text left-aligned)
#   (BOTTOM) 18D latent-dimension legend (boxed annotation; text left-aligned)
#
# LOCKED constraints preserved (NO geometry/camera/layout redesign beyond tasks):
#   - Geometry / hull construction / hull colors / barycenter logic unchanged.
#   - Vertex sign encoding (publication; circles only) is deterministic:
#       pos -> filled circle
#       neg -> circle-open
#     Implemented via TWO marker traces (pos + neg) to avoid per-point symbol collapse.
#   - Standalone HTML saving is deterministic under a folder rooted at getwd().
#   - mesh3d cleanup retained: no invalid `mode` attribute.
#
# ONLY CHANGES IN THIS VERSION (tasks requested):
#   (1) Increase axis-title ↔ tick-label separation using `title.standoff`.
#   (2) Increase font size uniformly for:
#         - Plotly trace legend,
#         - vertex/barycenter boxed annotation,
#         - 18D latent-dimension boxed annotation.
#   (3) Increase vertical separation between the TWO boxed annotations
#       (prevent touching/overlap during first paint).
#   (4) Stabilize first-open layout in saved standalone HTML via a post-load
#       Plotly resize + relayout hook (mitigates first-render annotation collisions
#       observed in some browsers).
###############################################################################

library(plotly)
library(dplyr)
library(htmlwidgets)
library(tools)

plot_circuitry_polytope_custom <- function(poly_obj,
                                           show_hull = TRUE,
                                           show_vertices = TRUE,
                                           show_barycenters = TRUE,
                                           title = NULL,
                                           outfile_html = NULL,
                                           outdir_html = "Figure_3_HTML") {
  
  verts_all <- poly_obj$vertices_all
  verts_sig <- poly_obj$vertices_sig
  verts_int <- poly_obj$vertices_int
  hull_sig  <- poly_obj$hull_sig
  hull_int  <- poly_obj$hull_int
  bary_sig  <- poly_obj$bary_sig_3d
  bary_int  <- poly_obj$bary_int_3d
  
  # -----------------------------
  # 18D dimension labels — authoritative vector
  # -----------------------------
  latent_dim_names <- poly_obj$latent_dim_names
  if (is.null(latent_dim_names) || length(latent_dim_names) != 18) {
    latent_dim_names <- c(
      "rho","rho_strength",
      "TN_dir","TN_strength",
      "OS_dir","OS_strength","OS_lr_chisq",
      "DSS_dir","DSS_strength","DSS_lr_chisq",
      "DFI_dir","DFI_strength","DFI_lr_chisq",
      "PFI_dir","PFI_strength","PFI_lr_chisq",
      "TME_score","Immune_dir"
    )
  }
  
  if (is.null(title)) {
    title <- paste0("Circuitry polytope (18D → 3D): ", poly_obj$circuitry_id)
  }
  
  p <- plot_ly()
  
  # -----------------------------
  # 1 — Vertices (LOCKED, circles-only sign encoding)
  # -----------------------------
  if (show_vertices && nrow(verts_all) > 0L) {
    
    stopifnot(all(verts_all$sign %in% c("pos","neg")))
    
    v_pos <- dplyr::filter(verts_all, sign == "pos")
    v_neg <- dplyr::filter(verts_all, sign == "neg")
    
    if (nrow(v_pos) > 0L) {
      p <- p %>%
        add_markers(
          data = v_pos,
          x = ~PC1, y = ~PC2, z = ~PC3,
          marker = list(size = 5, symbol = "circle", line = list(width = 1)),
          color  = ~side,
          colors = c("sig" = "darkblue", "int" = "darkorange"),
          text   = ~paste("Side:", side,
                          "<br>Dimension:", dim_name,
                          "<br>Sign:", sign),
          hoverinfo = "text",
          name = "Vertices (+)"
        )
    }
    
    if (nrow(v_neg) > 0L) {
      p <- p %>%
        add_markers(
          data = v_neg,
          x = ~PC1, y = ~PC2, z = ~PC3,
          marker = list(size = 5, symbol = "circle-open", line = list(width = 1)),
          color  = ~side,
          colors = c("sig" = "darkblue", "int" = "darkorange"),
          text   = ~paste("Side:", side,
                          "<br>Dimension:", dim_name,
                          "<br>Sign:", sign),
          hoverinfo = "text",
          name = "Vertices (−)"
        )
    }
  }
  
  # -----------------------------
  # 2 — Signature hull (baseline)
  # -----------------------------
  if (show_hull && !is.null(hull_sig$faces)) {
    
    simplices <- hull_sig$faces
    coords <- as.matrix(verts_sig[, c("PC1", "PC2", "PC3")])
    
    tri <- as.vector(t(simplices)) - 1L
    n_tri <- nrow(simplices)
    
    p <- p %>% add_trace(
      type   = "mesh3d",
      x      = coords[,1],
      y      = coords[,2],
      z      = coords[,3],
      i      = tri[seq(1, 3*n_tri, 3)],
      j      = tri[seq(2, 3*n_tri, 3)],
      k      = tri[seq(3, 3*n_tri, 3)],
      opacity     = 0.35,
      vertexcolor = rep("#2CA9BC", nrow(coords)),
      flatshading = TRUE,
      name        = "Signature hull"
    )
  }
  
  # -----------------------------
  # 3 — Interaction hull (baseline)
  # -----------------------------
  if (show_hull && !is.null(hull_int$faces)) {
    
    simplices <- hull_int$faces
    coords <- as.matrix(verts_int[, c("PC1", "PC2", "PC3")])
    
    tri <- as.vector(t(simplices)) - 1L
    n_tri <- nrow(simplices)
    
    p <- p %>% add_trace(
      type   = "mesh3d",
      x      = coords[,1],
      y      = coords[,2],
      z      = coords[,3],
      i      = tri[seq(1, 3*n_tri, 3)],
      j      = tri[seq(2, 3*n_tri, 3)],
      k      = tri[seq(3, 3*n_tri, 3)],
      opacity     = 0.35,
      vertexcolor = rep("#F2A65A", nrow(coords)),
      flatshading = TRUE,
      name        = "Interaction hull"
    )
  }
  
  # -----------------------------
  # 4 — Barycenters (baseline)
  # -----------------------------
  if (show_barycenters) {
    p <- p %>%
      add_markers(
        x = bary_sig["PC1"], y = bary_sig["PC2"], z = bary_sig["PC3"],
        marker = list(size = 10, symbol = "diamond", color = "darkblue"),
        name = "Signature barycenter"
      ) %>%
      add_markers(
        x = bary_int["PC1"], y = bary_int["PC2"], z = bary_int["PC3"],
        marker = list(size = 10, symbol = "diamond", color = "darkorange"),
        name = "Interaction barycenter"
      )
  }
  
  # -----------------------------
  # 5 — RIGHT legend column + axis pane background
  #     (Legend + two boxed annotations are positioned in paper coordinates;
  #      they do not depend on camera or scene interaction.)
  # -----------------------------
  LEGEND_COL_X0 <- 0.72
  LEGEND_X      <- LEGEND_COL_X0 + 0.02
  
  # Trace legend anchor (top of right column)
  TOP_Y <- 0.95
  
  # Task (3): vertical separation for boxed annotations (first-open safe spacing)
  # - Middle box is yanchor="top": increasing y moves it upward.
  # - Bottom box is yanchor="bottom": decreasing y moves it downward.
  MID_Y <- 0.695
  BOT_Y <- -0.05
  
  dim_legend_text <- paste0(
    "<b>Latent dimensions (18D):</b><br>",
    paste0(seq_along(latent_dim_names), " = ", latent_dim_names, collapse = "<br>")
  )
  
  vertex_legend_text <- paste0(
    "<b>Vertex legend:</b><br>",
    "Sign: + = filled circle; − = open circle<br>",
    "Side: <span style='color:darkblue;'><b>sig</b></span> = darkblue; ",
    "<span style='color:darkorange;'><b>int</b></span> = darkorange<br>",
    "Barycenters: diamonds (sig/int)"
  )
  
  # Axis pane + grid styling (Plotly 3D axes only; pane background is opaque light grey)
  PANE_BG_GREY <- "rgb(242,242,242)"
  GRID_COLOR   <- "rgba(0,0,0,0.45)"
  GRID_WIDTH   <- 2
  AXISLINE_COL <- "rgba(0,0,0,0.55)"
  ZERO_COL     <- "rgba(0,0,0,0.60)"
  
  # Task (1): increase axis-title ↔ tick-label separation (title.standoff)
  axis_with_grey_pane <- function(axis_title,
                                  axis_title_size = 25,
                                  tick_label_size = 12,
                                  title_standoff  = 90) {
    list(
      title = list(
        text = axis_title,
        font = list(size = axis_title_size),
        standoff = title_standoff
      ),
      tickfont = list(size = tick_label_size),
      showgrid = TRUE,
      gridcolor = GRID_COLOR,
      gridwidth = GRID_WIDTH,
      showbackground = TRUE,
      backgroundcolor = PANE_BG_GREY,
      showline = TRUE,
      linecolor = AXISLINE_COL,
      linewidth = 2,
      zeroline = TRUE,
      zerolinecolor = ZERO_COL,
      zerolinewidth = 2
    )
  }
  
  # Task (2): uniform legend-font scaling (trace legend + both boxed annotations)
  LEGEND_FONT_SCALE <- 1.25
  LEGEND_FONT_BASE  <- 11
  legend_font_size_scaled <- as.integer(round(LEGEND_FONT_BASE * LEGEND_FONT_SCALE))
  
  p <- p %>% layout(
    title = list(text = title),
    
    # Scene viewport is restricted to the left region to reserve the right column
    scene = list(
      domain = list(x = c(0.00, LEGEND_COL_X0), y = c(0.01, 1.00)),
      xaxis = axis_with_grey_pane("PC1"),
      yaxis = axis_with_grey_pane("PC2"),
      zaxis = axis_with_grey_pane("PC3")
    ),
    
    # Trace legend (Plotly-native legend; boxed)
    showlegend = TRUE,
    legend = list(
      x = LEGEND_X,
      xref = "paper",
      xanchor = "left",
      y = TOP_Y,
      yref = "paper",
      yanchor = "top",
      orientation = "v",
      traceorder = "normal",
      itemsizing = "constant",
      font = list(size = legend_font_size_scaled),
      bgcolor = "rgba(255,255,255,0.92)",
      bordercolor = "black",
      borderwidth = 0.5
    ),
    
    # Boxed annotations (paper coordinates; both are explicitly bordered and left-aligned)
    annotations = list(
      list(
        x = LEGEND_X, y = MID_Y,
        xref = "paper", yref = "paper",
        xanchor = "left", yanchor = "top",
        showarrow = FALSE,
        align = "left",
        text = vertex_legend_text,
        font = list(size = legend_font_size_scaled),
        bordercolor = "black",
        borderwidth = 0.5,
        bgcolor = "rgba(255,255,255,0.92)"
      ),
      list(
        x = LEGEND_X, y = BOT_Y,
        xref = "paper", yref = "paper",
        xanchor = "left", yanchor = "bottom",
        showarrow = FALSE,
        align = "left",
        text = dim_legend_text,
        font = list(size = legend_font_size_scaled),
        bordercolor = "black",
        borderwidth = 0.5,
        bgcolor = "rgba(255,255,255,0.92)"
      )
    ),
    
    # Layout margins (outer whitespace around the full widget):
    # l/r = room for axes + right column; t = title clearance; b = bottom clearance
    margin = list(l = 60, r = 10, t = 40, b = 40)
  )
  
  # -----------------------------
  # 6 — Save HTML if requested (deterministic folder under wd + verify)
  #     Task (4): stabilize first-open layout in standalone HTML
  # -----------------------------
  if (!is.null(outfile_html)) {
    
    # Canvas size constants are retained for documentation/consistency only.
    # (No sizingPolicy wrapper is applied in this locked baseline.)
    FIG_W <- 2600
    FIG_H <- 1600
    
    wd_abs <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    outdir_abs <- file.path(wd_abs, outdir_html)
    
    if (!dir.exists(outdir_abs)) {
      dir.create(outdir_abs, recursive = TRUE, showWarnings = FALSE)
    }
    stopifnot(dir.exists(outdir_abs))
    
    target_html <- file.path(outdir_abs, outfile_html)
    if (tolower(tools::file_ext(target_html)) != "html") {
      target_html <- paste0(target_html, ".html")
    }
    
    # No sizingPolicy wrapper (locked baseline). Only a post-load redraw hook is applied.
    p_sized <- p
    
    # Standalone-HTML stabilization:
    # Force Plotly to finalize layout after DOM + fonts + WebGL sizing,
    # reducing the risk of initial annotation/legend collisions on first paint.
    p_sized <- htmlwidgets::onRender(
      p_sized,
      "
          function(el, x) {
            function relayoutOnce() {
              try {
                if (window.Plotly && el && el.id) {
                Plotly.Plots.resize(el);
               Plotly.relayout(el, {'scene.domain.y': [0.01, 1.00]});
                }
              } catch(e) {}
            }
    
            // Run once shortly after render, and again after a brief delay
            // to capture font loading + initial WebGL sizing.
            setTimeout(relayoutOnce, 50);
            setTimeout(relayoutOnce, 400);
    
            // Re-apply once on window load (browser timing edge cases).
            window.addEventListener('load', function(){ setTimeout(relayoutOnce, 50); });
          }
          "
    )
    
    htmlwidgets::saveWidget(
      widget        = p_sized,
      file          = target_html,
      selfcontained = TRUE
    )
    
    stopifnot(file.exists(target_html))
    message("OK — HTML saved: ",
            normalizePath(target_html, winslash = "/", mustWork = TRUE))
  }
  
  return(p)
}

###############################################################################
# USAGE EXAMPLE (MUST be outside the function)
###############################################################################

# setwd("D:/Pré-artigo 5-optosis model/Omica signatures the multidimensional concept")

poly_LGG <- build_circuitry_polytope(
  tensor, embedding,
  circuitry_id = "LGG-3207 / LGG-2427"
)

p <- plot_circuitry_polytope_custom(
  poly_LGG,
  outfile_html = "Final_LGG_3207_2427_custom.html",
  outdir_html  = "Figure_3_HTML"
)

print(p)

### 
### 
### 
### OTHER USAGES EXAMPLES
### 
### 
### 

### Other example based on the universla baseline reference function 
poly_PAAD <- build_circuitry_polytope(
  tensor, embedding,
  circuitry_id = "PAAD-3581 / PAAD-7300"
)

p <- plot_circuitry_polytope_custom(
  poly_PAAD,
  outfile_html = "Final_PAAD_3581_7300_custom.html",
  outdir_html  = "Figure_N_HTML"
)

print(p)

### Other example based on the universal baseline reference function 
poly_PAAD <- build_circuitry_polytope(
  tensor, embedding,
  circuitry_id = "PAAD-3903 / PAAD-7083"
)

p <- plot_circuitry_polytope_custom(
  poly_PAAD,
  outfile_html = "Final_PAAD_3903_7083_custom.html",
  outdir_html  = "Figure_N_HTML"
)

print(p)

### Other example based on the universal baseline reference function 
poly_PAAD <- build_circuitry_polytope(
  tensor, embedding,
  circuitry_id = "PAAD-10164 / PAAD-2159"
)

p <- plot_circuitry_polytope_custom(
  poly_PAAD,
  outfile_html = "Final_PAAD_10164_2179_custom.html",
  outdir_html  = "Figure_N_HTML"
)

print(p)


###############################################################################
# Figure 3 — RASTERIZATION PIPELINE (HTML → headless Chrome PNG → 600 DPI TIFF)
#
# Tuning goal (current):
#   - Text larger than the first “tiny-text” attempt
#   - Avoid boxed-legend overlap (superimposition)
#
# Strategy:
#   - Moderate CSS viewport (not too small)
#   - High DSF for pixel resolution
#
# Current tuned settings:
#   CSS viewport 1100×850 with DSF=5 → PNG ≈ 5500×4250 pixels
###############################################################################

suppressPackageStartupMessages({
  library(processx)
  library(magick)
})

############################
# USER-ADJUSTABLE SETTINGS #
############################

HTML_IN <- "Figure_3_HTML/Final_LGG_3207_2427_custom.html"

PNG_OUT <- "Figure_3_HTML/Final_LGG_3207_2427_highres.png"
TIF_OUT <- "Figure_3_HTML/Final_LGG_3207_2427_600dpi.tiff"

# --- Layout (controls relative text size and overlap pressure) ----------------
CSS_VIEWPORT_WIDTH  <- 1100
CSS_VIEWPORT_HEIGHT <- 850

# --- Resolution (controls pixels, not relative sizing) ------------------------
DEVICE_SCALE_FACTOR <- 5
# ------------------------------------------------------------------------------

WAIT_AFTER_LOAD_MS <- 9000
CHROME_TIMEOUT_SEC <- 180

TARGET_DPI <- 600
CHROME_LOG <- "chrome_headless_log.txt"

###########################
# RESOLVE ABSOLUTE PATHS  #
###########################

html_in_abs <- normalizePath(HTML_IN, winslash = "/", mustWork = TRUE)
png_out_abs <- normalizePath(PNG_OUT, winslash = "/", mustWork = FALSE)
tif_out_abs <- normalizePath(TIF_OUT, winslash = "/", mustWork = FALSE)
chrome_log_abs <- normalizePath(CHROME_LOG, winslash = "/", mustWork = FALSE)

message("HTML  : ", html_in_abs)
message("PNG   : ", png_out_abs)
message("TIFF  : ", tif_out_abs)
message("LOG   : ", chrome_log_abs)

dir.create(dirname(png_out_abs), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(tif_out_abs), recursive = TRUE, showWarnings = FALSE)

#################################
# LOCATE CHROME / CHROMIUM (WIN) #
#################################

find_chrome <- function() {
  candidates <- c(
    Sys.getenv("CHROME_PATH", unset = ""),
    "C:/Program Files/Google/Chrome/Application/chrome.exe",
    "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
    "C:/Program Files/Chromium/Application/chrome.exe",
    "C:/Program Files (x86)/Chromium/Application/chrome.exe"
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

chrome_exe <- find_chrome()
if (is.na(chrome_exe)) {
  stop(
    "Chrome/Chromium not found.\n",
    "Install Google Chrome or set CHROME_PATH to chrome.exe."
  )
}
message("Chrome detected: ", chrome_exe)

###############################################################################
# STAGING: ASCII-ONLY TEMP WORKSPACE (prevents Windows file:/// Unicode issues)
###############################################################################

stage_root <- file.path(
  tempdir(),
  paste0("fig3_raster_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)

profile_dir <- file.path(stage_root, "chrome_profile")
dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)

html_stage_orig <- file.path(stage_root, "figure3.html")
ok_copy <- file.copy(from = html_in_abs, to = html_stage_orig, overwrite = TRUE)
if (!isTRUE(ok_copy)) stop("Failed to stage HTML into temp directory: ", html_stage_orig)

png_stage <- file.path(stage_root, "figure3.png")
log_stage <- file.path(stage_root, "chrome_log.txt")

if (file.exists(png_stage)) file.remove(png_stage)
if (file.exists(log_stage)) file.remove(log_stage)

# Wrapper to force full-viewport sizing
html_stage_wrap <- file.path(stage_root, "wrapper.html")

orig_url <- paste0(
  "file:///",
  gsub("\\\\", "/", normalizePath(html_stage_orig, winslash = "/", mustWork = TRUE))
)

wrapper_txt <- sprintf(
  '<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <style>
    html, body { margin:0; padding:0; width:100%%; height:100%%; overflow:hidden; background:#ffffff; }
    iframe { border:0; width:100%%; height:100%%; display:block; }
  </style>
</head>
<body>
  <iframe src="%s"></iframe>
</body>
</html>',
orig_url
)

writeLines(wrapper_txt, html_stage_wrap)

file_url <- paste0(
  "file:///",
  gsub("\\\\", "/", normalizePath(html_stage_wrap, winslash = "/", mustWork = TRUE))
)

message("STAGE DIR : ", stage_root)
message("WRAPPER  : ", html_stage_wrap)
message("STAGE URL : ", file_url)
message(
  "CSS viewport: ", CSS_VIEWPORT_WIDTH, " x ", CSS_VIEWPORT_HEIGHT,
  " ; DSF=", DEVICE_SCALE_FACTOR,
  " ; expected PNG ≈ ",
  CSS_VIEWPORT_WIDTH * DEVICE_SCALE_FACTOR, " x ",
  CSS_VIEWPORT_HEIGHT * DEVICE_SCALE_FACTOR, " px"
)

############################################
# HEADLESS CHROME SCREENSHOT               #
############################################

chrome_args <- c(
  "--headless=new",
  "--disable-gpu",
  "--hide-scrollbars",
  
  paste0("--user-data-dir=", profile_dir),
  
  "--no-first-run",
  "--no-default-browser-check",
  "--disable-background-networking",
  "--disable-sync",
  "--metrics-recording-only",
  "--disable-features=Translate,MediaRouter",
  
  "--allow-file-access-from-files",
  
  sprintf("--window-size=%d,%d", CSS_VIEWPORT_WIDTH, CSS_VIEWPORT_HEIGHT),
  sprintf("--force-device-scale-factor=%d", DEVICE_SCALE_FACTOR),
  
  sprintf("--virtual-time-budget=%d", WAIT_AFTER_LOAD_MS),
  "--run-all-compositor-stages-before-draw",
  
  paste0("--screenshot=", png_stage),
  
  file_url
)

p <- processx::process$new(
  command = chrome_exe,
  args    = chrome_args,
  stdout  = log_stage,
  stderr  = log_stage
)

t0 <- Sys.time()
repeat {
  Sys.sleep(0.2)
  
  if (file.exists(png_stage)) {
    sz <- file.info(png_stage)$size
    if (!is.na(sz) && sz > 50e3) break
  }
  
  if (!p$is_alive()) break
  
  if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > CHROME_TIMEOUT_SEC) {
    p$kill()
    break
  }
}

if (p$is_alive()) p$kill()

if (!file.exists(png_stage)) {
  cat("\nHeadless Chrome did not produce the staged PNG within timeout.\n")
  cat("Staged log file: ", log_stage, "\n", sep = "")
  cat("\n---- Chrome staged log (last ~4000 chars) ----\n")
  if (file.exists(log_stage)) {
    txt <- paste(readLines(log_stage, warn = FALSE), collapse = "\n")
    cat(substr(txt, max(1, nchar(txt) - 4000), nchar(txt)), "\n")
  } else {
    cat("(No staged log file created)\n")
  }
  if (file.exists(log_stage)) file.copy(log_stage, chrome_log_abs, overwrite = TRUE)
  stop("Headless capture failed/stalled. Inspect staged log and/or chrome_headless_log.txt.")
}

message("Staged PNG created: ", png_stage, " (", file.info(png_stage)$size, " bytes)")

file.copy(from = png_stage, to = png_out_abs, overwrite = TRUE)
if (file.exists(log_stage)) file.copy(from = log_stage, to = chrome_log_abs, overwrite = TRUE)
message("Final PNG copied: ", png_out_abs)

info_png <- magick::image_info(magick::image_read(png_out_abs))
message("PNG dimensions: ", info_png$width, " x ", info_png$height, " px")

message("Final PNG copied: ", png_out_abs)

# ------------------------------------------------------------------
# OPTIONAL: Embed 600 dpi metadata into PNG (pixels unchanged)
# NOTE:
# - This does NOT resample or change resolution
# - Many viewers (e.g. Windows Explorer) still show 96 dpi
# - Harmless and sometimes useful for downstream tools
# ------------------------------------------------------------------
img_png <- magick::image_read(png_out_abs)

magick::image_write(
  img_png,
  path   = png_out_abs,
  format = "png",
  density = "600x600"
)


############################################
# PNG -> 600 DPI TIFF (NO RESAMPLING)      #
############################################

img <- magick::image_read(png_out_abs)
density_str <- paste0(TARGET_DPI, "x", TARGET_DPI)

magick::image_write(
  img,
  path = tif_out_abs,
  format = "tiff",
  compression = "lzw",
  density = density_str
)

message("TIFF created: ", tif_out_abs, " (density metadata = ", density_str, ")")

library(magick)

# Check pixel dimensions (this is the real resolution)
info_png <- image_info(image_read("Figure_3_HTML/Final_LGG_3207_2427_highres.png"))
print(info_png[, c("width","height")])

# Check TIFF density tag (this is the 600 dpi metadata)
info_tif <- image_info(image_read("Figure_3_HTML/Final_LGG_3207_2427_600dpi.tiff"))
print(info_tif)
  
