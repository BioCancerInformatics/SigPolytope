# ============================================================
# R/build_tensor_semantic.R
# ============================================================
# Global semantic encoder -> tensor (features_sig/features_int)
# Works for any data.frame as long as columns follow recognizable patterns.
# Does NOT assume any specific dataset file.
# ============================================================

#' @keywords internal
.default_dimension_schema <- function() {
  list(
    # -----------------------------
    # Tumor vs Normal
    # -----------------------------
    TN_dir = list(
      kind = "dir_categorical",
      pattern = "Tumor_vs_normal_type",  # sem _sig/_int ainda
      encoder = function(x) {
        z <- tolower(trimws(as.character(x)))
        dplyr::case_when(
          z %in% c("overexpression", "over", "up", "upregulated") ~  1,
          z %in% c("underexpression", "under", "down", "downregulated") ~ -1,
          z %in% c("unchanged", "ns", "no data", "nodata", "na", "n/a", "") ~ 0,
          is.na(z) ~ 0,
          TRUE ~ 0
        )
      }
    ),

    TN_strength = list(
      kind = "numeric_strength",
      pattern = "Tumor_vs_normal_strength",
      encoder = function(x) suppressWarnings(as.numeric(x))
    ),

    # -----------------------------
    # Survival endpoints (generic template, expanded per endpoint)
    # -----------------------------
    SURV_dir = list(
      kind = "dir_survival",
      pattern = "Cox_(OS|DSS|DFI|PFI)_type",
      encoder = function(x) {
        z <- tolower(trimws(as.character(x)))
        dplyr::case_when(
          z %in% c("risk", "hazard", "hr>1", "highrisk", "high risk") ~  1,
          z %in% c("protective", "hr<1", "lowrisk", "low risk") ~ -1,
          z %in% c("ns", "no data", "nodata", "na", "n/a", "") ~ 0,
          is.na(z) ~ 0,
          TRUE ~ 0
        )
      }
    ),

    SURV_strength = list(
      kind = "numeric_strength",
      pattern = "Cox_(OS|DSS|DFI|PFI)_strength",
      encoder = function(x) suppressWarnings(as.numeric(x))
    ),

    SURV_chi2 = list(
      kind = "numeric_chi2",
      # aceitar chisq ou chi2 (variações comuns)
      pattern = "Cox_(OS|DSS|DFI|PFI)_(chisq|chi2)",
      encoder = function(x) suppressWarnings(as.numeric(x))
    )
  )
}

#' @keywords internal
.find_first_matching_col <- function(cols, pattern, suffix) {
  # pattern pode ser regex; o sufixo define o lado (_sig/_int)
  rx <- paste0("^", pattern, suffix, "$")
  hit <- grep(rx, cols, value = TRUE)
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

#' @keywords internal
.check_required_groups <- function(data, side_suffix, require = c("TN", "SURV")) {
  cn <- colnames(data)

  has_any <- function(rx) any(grepl(rx, cn))

  if ("TN" %in% require) {
    ok_sig <- has_any(paste0("^Tumor_vs_normal_(type|strength)", side_suffix[["sig"]], "$"))
    ok_int <- has_any(paste0("^Tumor_vs_normal_(type|strength)", side_suffix[["int"]], "$"))
    if (!(ok_sig && ok_int)) {
      warning("Grupo TN incompleto em sig/int (pacote global: prosseguindo).")
    }
  }

  if ("SURV" %in% require) {
    ok_sig <- has_any(paste0("^Cox_[A-Za-z0-9]+_(type|strength|chisq|chi2)", side_suffix[["sig"]], "$"))
    ok_int <- has_any(paste0("^Cox_[A-Za-z0-9]+_(type|strength|chisq|chi2)", side_suffix[["int"]], "$"))
    if (!(ok_sig && ok_int)) {
      warning("Grupo SURV incompleto em sig/int (pacote global: prosseguindo).")
    }
  }

  invisible(TRUE)
}

#' @keywords internal
.expand_schema_by_data <- function(schema, data, side_suffix,
                                   mode = c("discovered", "canonical"),
                                   canonical_endpoints = c("OS", "DSS", "DFI", "PFI")) {
  stopifnot(is.list(schema), is.data.frame(data))
  mode <- match.arg(mode)

  cn <- colnames(data)

  # Detecta endpoints presentes em ambos os lados (sig e int)
  rx_sig <- paste0("^Cox_([A-Za-z0-9]+)_(type|strength|chisq|chi2)", side_suffix[["sig"]], "$")
  rx_int <- paste0("^Cox_([A-Za-z0-9]+)_(type|strength|chisq|chi2)", side_suffix[["int"]], "$")

  get_eps <- function(rx) {
    m <- regexec(rx, cn)
    hits <- regmatches(cn, m)
    unique(vapply(hits, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)))
  }

  eps <- unique(c(get_eps(rx_sig), get_eps(rx_int)))
  eps <- eps[!is.na(eps)]

  endpoints <- if (mode == "canonical") {
    unique(c(canonical_endpoints, eps))
  } else {
    eps
  }

  out <- list()

  for (nm in names(schema)) {
    spec <- schema[[nm]]

    # Expande apenas specs com o bloco (OS|DSS|DFI|PFI)
    if (grepl("\\(OS\\|DSS\\|DFI\\|PFI\\)", spec$pattern)) {
      for (ep in endpoints) {
        nm2  <- sub("^SURV_", paste0(ep, "_"), nm)  # SURV_dir -> OS_dir, etc.
        pat2 <- gsub("\\(OS\\|DSS\\|DFI\\|PFI\\)", ep, spec$pattern)

        out[[nm2]] <- spec
        out[[nm2]]$pattern  <- pat2
        out[[nm2]]$endpoint <- ep
      }
    } else {
      out[[nm]] <- spec
    }
  }

  # ------------------------------------------------------------
  # ✅ FIX: Ordenação canônica SEM descartar blocos extras
  # ------------------------------------------------------------
  # (1) uma ordem "preferida" (se existir):
  ord <- c(
    grep("^CORR_",  names(out), value = TRUE),
    grep("^TN_",    names(out), value = TRUE),
    unlist(lapply(endpoints, function(ep) grep(paste0("^", ep, "_"), names(out), value = TRUE))),
    grep("^MICRO_", names(out), value = TRUE)
  )

  # (2) anexa tudo que sobrou (não perde dimensões!)
  rest <- setdiff(names(out), ord)

  out[unique(c(ord, rest))]
}


#' @keywords internal
.encode_side <- function(data, schema, suffix, mode = c("discovered", "canonical")) {
  mode <- match.arg(mode)

  feats <- list()
  cols_used <- list()
  n <- nrow(data)

  for (dim_name in names(schema)) {
    spec <- schema[[dim_name]]
    col  <- .find_first_matching_col(colnames(data), spec$pattern, suffix)

    if (is.null(col)) {
      if (mode == "canonical") {
        feats[[dim_name]] <- rep(0, n)
        cols_used[[dim_name]] <- NA_character_
      }
      next
    }

    v <- data[[col]]
    x <- spec$encoder(v)

    # blindagem: garante numeric e estabilidade
    x <- suppressWarnings(as.numeric(x))
    x[!is.finite(x)] <- 0

    feats[[dim_name]] <- x
    cols_used[[dim_name]] <- col
  }

  list(
    features = as.data.frame(feats, check.names = FALSE),
    columns  = unlist(cols_used, use.names = TRUE)
  )
}

#' @keywords internal
.apply_dir_mask_zeroing <- function(features_df) {
  # Regra global:
  # - se TN_dir == 0 -> TN_strength = 0
  # - se <EP>_dir == 0 -> <EP>_strength = 0 e <EP>_chi2 = 0
  if (is.null(features_df) || ncol(features_df) == 0) return(features_df)

  cn <- colnames(features_df)

  # TN
  if (all(c("TN_dir", "TN_strength") %in% cn)) {
    idx <- which(features_df[["TN_dir"]] == 0)
    if (length(idx) > 0) features_df[["TN_strength"]][idx] <- 0
  }

  # endpoints = prefixos antes de "_dir"
  dir_cols <- grep("^[A-Za-z0-9]+_dir$", cn, value = TRUE)
  endpoints <- sub("_dir$", "", dir_cols)

  # evita reprocessar TN (já feito acima)
  endpoints <- setdiff(endpoints, "TN")

  for (ep in endpoints) {
    dir_name <- paste0(ep, "_dir")
    strength <- paste0(ep, "_strength")
    chi2     <- paste0(ep, "_chi2")

    if (!dir_name %in% cn) next

    idx <- which(features_df[[dir_name]] == 0)
    if (length(idx) == 0) next

    if (strength %in% cn) features_df[[strength]][idx] <- 0
    if (chi2     %in% cn) features_df[[chi2]][idx] <- 0
  }

  features_df
}

#' Build geometric tensor from semantic columns (global encoder)
#'
#' Constrói automaticamente `features_sig` e `features_int` a partir de um
#' data.frame com colunas semânticas (ex.: Cox_OS_type_sig, Tumor_vs_normal_strength_int, etc.).
#' O pacote permanece global: as dimensões são descobertas via padrões de nomes.
#'
#' @param data data.frame/tibble.
#' @param id_col Nome da coluna de ID (default "Circuitries_id").
#' @param side_suffix Named character vector com sufixos dos lados.
#' @param schema Schema global (lista) com patterns e encoders.
#' @param min_dims Mínimo de dimensões comuns exigidas.
#' @param mode "discovered" usa apenas dimensões encontradas; "canonical" mantém dimensões-alvo,
#'   preenchendo ausentes com 0.
#' @param canonical_endpoints Endpoints padrão usados no modo "canonical".
#'
#' @return Lista com: features_sig, features_int, meta, latent_dims.
#' @export
build_geometric_tensor <- function(
    data,
    id_col = "Circuitries_id",
    side_suffix = c(sig = "_sig", int = "_int"),
    schema = .default_dimension_schema(),
    min_dims = 3,
    mode = c("discovered", "canonical"),
    canonical_endpoints = c("OS", "DSS", "DFI", "PFI")
) {
  stopifnot(is.data.frame(data))
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Precisa de dplyr.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Precisa de tibble.")

  mode <- match.arg(mode)

  # garante side_suffix nomeado (evita bug silencioso)
  if (is.null(names(side_suffix)) || !all(c("sig", "int") %in% names(side_suffix))) {
    stop("side_suffix precisa ser nomeado: c(sig='_sig', int='_int').")
  }

  .check_required_groups(data, side_suffix)

  meta <- if (id_col %in% colnames(data)) data[, id_col, drop = FALSE] else data.frame()

  schema2 <- .expand_schema_by_data(
    schema = schema,
    data = data,
    side_suffix = side_suffix,
    mode = mode,
    canonical_endpoints = canonical_endpoints
  )

  out_sig <- .encode_side(data, schema2, suffix = side_suffix[["sig"]], mode = mode)
  out_int <- .encode_side(data, schema2, suffix = side_suffix[["int"]], mode = mode)

  # ✅ aplica regra global de coerência: dir==0 zera strength/chi2
  out_sig$features <- .apply_dir_mask_zeroing(out_sig$features)
  out_int$features <- .apply_dir_mask_zeroing(out_int$features)

  dims_keep <- intersect(names(out_sig$features), names(out_int$features))

  # canonical: mantém a ordem canônica do schema2
  if (mode == "canonical") {
    dims_keep <- names(schema2)[names(schema2) %in% dims_keep]
  }

  if (length(dims_keep) < min_dims) {
    stop("Poucas dimensões comuns entre sig/int: ", length(dims_keep), ".")
  }

  X_sig <- as.matrix(out_sig$features[, dims_keep, drop = FALSE])
  X_int <- as.matrix(out_int$features[, dims_keep, drop = FALSE])

  list(
    features_sig = X_sig,
    features_int = X_int,
    meta = tibble::as_tibble(meta),
    latent_dims = list(
      canonical   = dims_keep,
      columns_sig = out_sig$columns[dims_keep],
      columns_int = out_int$columns[dims_keep],
      mode        = mode
    )
  )
}
