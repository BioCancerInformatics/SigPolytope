# ============================================================
# RUN GeoCircuitry — 18D PAPER-ALIGNED LATENT SPACE (FIXED)
# - builds 18D features for sig/int with IDENTICAL colnames
# - avoids pairs API
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(SigPolytope)
})
# setwd("G:/Projetos_V2/SigPolytope")
# devtools::load_all()

# library(SigPolytope)


# setwd("G:/Projetos_V2/SigPolytope")  # local do pacote
# devtools::load_all()

setwd("G:/Pacote geometria")          # local do TSV

# -----------------------------
# Helpers
# -----------------------------
coerce_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(x))
}

neg_log10_p <- function(p, eps = 1e-12) {
  p_num <- coerce_numeric_safe(p)
  out <- rep(NA_real_, length(p_num))
  ok <- is.finite(p_num) & p_num > 0
  out[ok] <- -log10(p_num[ok] + eps)
  out
}

sanitize_matrix <- function(X, fill = 0) {
  X <- as.matrix(X)
  X[!is.finite(X)] <- fill
  X
}

encode_tvn_dir <- function(x) {
  z <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    z %in% c("overexpression","over","up","upregulated") ~  1,
    z %in% c("underexpression","under","down","downregulated") ~ -1,
    z %in% c("unchanged","ns","no data","nodata","na","n/a","") ~ 0,
    is.na(z) ~ 0,
    TRUE ~ 0
  )
}

encode_surv_dir <- function(x) {
  z <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    z %in% c("risk","hazard","hr>1","highrisk","high risk") ~  1,
    z %in% c("protective","hr<1","lowrisk","low risk") ~ -1,
    z %in% c("ns","no data","nodata","na","n/a","") ~ 0,
    is.na(z) ~ 0,
    TRUE ~ 0
  )
}

# NOTE: adjust to your real labels if needed
encode_immune_dir <- function(x) {
  z <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    z %in% c("hot","inflamed") ~  1,
    z %in% c("cold","desert")  ~ -1,
    z %in% c("suppressive","excluded") ~ 0,
    z %in% c("ns","no data","nodata","na","n/a","") ~ 0,
    is.na(z) ~ 0,
    TRUE ~ 0
  )
}

# -----------------------------
# 1) Load data
# -----------------------------
df <- read_tsv("Regulatory_circuitries.tsv", show_col_types = FALSE)

stopifnot("Circuitries_id" %in% colnames(df))

# -----------------------------
# 2) Required cols check (raw columns)
# -----------------------------
required_cols <- c(
  "Correlation_rho_sig", "Correlation_p.adj_sig",
  "Correlation_rho_int", "Correlation_p.adj_int",
  "Tumor_vs_normal_sig", "Tumor_vs_normal_p.adj_sig",
  "Tumor_vs_normal_int", "Tumor_vs_normal_p.adj_int",
  "Cox_OS_type_sig",  "Cox_OS_p.value_sig",  "OS_log_rank_chisq_sig",
  "Cox_OS_type_int",  "Cox_OS_p.value_int",  "OS_log_rank_chisq_int",
  "Cox_DSS_type_sig", "Cox_DSS_p.value_sig", "DSS_log_rank_chisq_sig",
  "Cox_DSS_type_int", "Cox_DSS_p.value_int", "DSS_log_rank_chisq_int",
  "Cox_DFI_type_sig", "Cox_DFI_p.value_sig", "DFI_log_rank_chisq_sig",
  "Cox_DFI_type_int", "Cox_DFI_p.value_int", "DFI_log_rank_chisq_int",
  "Cox_PFI_type_sig", "Cox_PFI_p.value_sig", "PFI_log_rank_chisq_sig",
  "Cox_PFI_type_int", "Cox_PFI_p.value_int", "PFI_log_rank_chisq_int",
  "Microenvironment_score_sig", "Microenvironment_score_int",
  "Immune_classification_sig",  "Immune_classification_int"
)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Colunas ausentes no TSV:\n- ", paste(missing_cols, collapse = "\n- "))
}

# -----------------------------
# 3) Build raw 18D components with _sig/_int suffix (still OK)
# -----------------------------
df2 <- df %>%
  mutate(
    # correlation
    rho_sig = coerce_numeric_safe(Correlation_rho_sig),
    rho_int = coerce_numeric_safe(Correlation_rho_int),
    rho_strength_sig = neg_log10_p(Correlation_p.adj_sig),
    rho_strength_int = neg_log10_p(Correlation_p.adj_int),
    
    # TN
    TN_dir_sig = encode_tvn_dir(Tumor_vs_normal_sig),
    TN_dir_int = encode_tvn_dir(Tumor_vs_normal_int),
    TN_strength_sig = neg_log10_p(Tumor_vs_normal_p.adj_sig),
    TN_strength_int = neg_log10_p(Tumor_vs_normal_p.adj_int),
    
    # OS
    OS_dir_sig      = encode_surv_dir(Cox_OS_type_sig),
    OS_dir_int      = encode_surv_dir(Cox_OS_type_int),
    OS_strength_sig = neg_log10_p(Cox_OS_p.value_sig),
    OS_strength_int = neg_log10_p(Cox_OS_p.value_int),
    OS_lr_chisq_sig = coerce_numeric_safe(OS_log_rank_chisq_sig),
    OS_lr_chisq_int = coerce_numeric_safe(OS_log_rank_chisq_int),
    
    # DSS
    DSS_dir_sig      = encode_surv_dir(Cox_DSS_type_sig),
    DSS_dir_int      = encode_surv_dir(Cox_DSS_type_int),
    DSS_strength_sig = neg_log10_p(Cox_DSS_p.value_sig),
    DSS_strength_int = neg_log10_p(Cox_DSS_p.value_int),
    DSS_lr_chisq_sig = coerce_numeric_safe(DSS_log_rank_chisq_sig),
    DSS_lr_chisq_int = coerce_numeric_safe(DSS_log_rank_chisq_int),
    
    # DFI
    DFI_dir_sig      = encode_surv_dir(Cox_DFI_type_sig),
    DFI_dir_int      = encode_surv_dir(Cox_DFI_type_int),
    DFI_strength_sig = neg_log10_p(Cox_DFI_p.value_sig),
    DFI_strength_int = neg_log10_p(Cox_DFI_p.value_int),
    DFI_lr_chisq_sig = coerce_numeric_safe(DFI_log_rank_chisq_sig),
    DFI_lr_chisq_int = coerce_numeric_safe(DFI_log_rank_chisq_int),
    
    # PFI
    PFI_dir_sig      = encode_surv_dir(Cox_PFI_type_sig),
    PFI_dir_int      = encode_surv_dir(Cox_PFI_type_int),
    PFI_strength_sig = neg_log10_p(Cox_PFI_p.value_sig),
    PFI_strength_int = neg_log10_p(Cox_PFI_p.value_int),
    PFI_lr_chisq_sig = coerce_numeric_safe(PFI_log_rank_chisq_sig),
    PFI_lr_chisq_int = coerce_numeric_safe(PFI_log_rank_chisq_int),
    
    # TME + immune
    TME_score_sig = coerce_numeric_safe(Microenvironment_score_sig),
    TME_score_int = coerce_numeric_safe(Microenvironment_score_int),
    
    Immune_dir_sig = encode_immune_dir(Immune_classification_sig),
    Immune_dir_int = encode_immune_dir(Immune_classification_int)
  )

# -----------------------------
# 4) Build two feature matrices with IDENTICAL column names (canonical 18D)
# -----------------------------
cols_18 <- c(
  "rho",
  "rho_strength",
  "TN_dir",
  "TN_strength",
  "OS_dir",
  "OS_strength",
  "OS_lr_chisq",
  "DSS_dir",
  "DSS_strength",
  "DSS_lr_chisq",
  "DFI_dir",
  "DFI_strength",
  "DFI_lr_chisq",
  "PFI_dir",
  "PFI_strength",
  "PFI_lr_chisq",
  "TME_score",
  "Immune_dir"
)

X_sig <- df2 %>%
  transmute(
    rho          = rho_sig,
    rho_strength = rho_strength_sig,
    TN_dir       = TN_dir_sig,
    TN_strength  = TN_strength_sig,
    OS_dir       = OS_dir_sig,
    OS_strength  = OS_strength_sig,
    OS_lr_chisq  = OS_lr_chisq_sig,
    DSS_dir      = DSS_dir_sig,
    DSS_strength = DSS_strength_sig,
    DSS_lr_chisq = DSS_lr_chisq_sig,
    DFI_dir      = DFI_dir_sig,
    DFI_strength = DFI_strength_sig,
    DFI_lr_chisq = DFI_lr_chisq_sig,
    PFI_dir      = PFI_dir_sig,
    PFI_strength = PFI_strength_sig,
    PFI_lr_chisq = PFI_lr_chisq_sig,
    TME_score    = TME_score_sig,
    Immune_dir   = Immune_dir_sig
  ) %>%
  as.matrix()

X_int <- df2 %>%
  transmute(
    rho          = rho_int,
    rho_strength = rho_strength_int,
    TN_dir       = TN_dir_int,
    TN_strength  = TN_strength_int,
    OS_dir       = OS_dir_int,
    OS_strength  = OS_strength_int,
    OS_lr_chisq  = OS_lr_chisq_int,
    DSS_dir      = DSS_dir_int,
    DSS_strength = DSS_strength_int,
    DSS_lr_chisq = DSS_lr_chisq_int,
    DFI_dir      = DFI_dir_int,
    DFI_strength = DFI_strength_int,
    DFI_lr_chisq = DFI_lr_chisq_int,
    PFI_dir      = PFI_dir_int,
    PFI_strength = PFI_strength_int,
    PFI_lr_chisq = PFI_lr_chisq_int,
    TME_score    = TME_score_int,
    Immune_dir   = Immune_dir_int
  ) %>%
  as.matrix()

colnames(X_sig) <- cols_18
colnames(X_int) <- cols_18

# sanitize
X_sig <- sanitize_matrix(X_sig, fill = 0)
X_int <- sanitize_matrix(X_int, fill = 0)

# -----------------------------
# 5) Build sig_geom/reg_geom with shared colnames
# -----------------------------
meta_cols <- setdiff(colnames(df2), c(
  # all computed raw columns + raw required cols
  colnames(df2 %>% select(matches("^(rho|rho_strength|TN_|OS_|DSS_|DFI_|PFI_|TME_|Immune_)")))
))

# Use a plain data.frame for meta to avoid tibble rownames warning
meta_df <- df2[, intersect(c("Circuitries_id", meta_cols), colnames(df2)), drop = FALSE]
meta_df <- as.data.frame(meta_df)
rownames(meta_df) <- df2$Circuitries_id

sig_geom <- structure(
  list(
    features  = X_sig,
    meta      = tibble::as_tibble(meta_df),
    dim_names = colnames(X_sig)
  ),
  class = "signature_geometry"
)
reg_geom <- structure(
  list(
    features  = X_int,
    meta      = tibble::as_tibble(meta_df),
    dim_names = colnames(X_int)
  ),
  class = c("regulator_geometry", "signature_geometry")
)

rownames(sig_geom$features) <- df2$Circuitries_id
rownames(reg_geom$features) <- df2$Circuitries_id

stopifnot(identical(colnames(sig_geom$features), colnames(reg_geom$features)))
stopifnot(ncol(sig_geom$features) == 18L)

# -----------------------------
# 6) Convergence + plot one circuitry
# -----------------------------
conv <- compute_circuitry_convergence(
  sig_geom, reg_geom,
  n_components = 3,
  parallel = FALSE,
  strict_dims = TRUE
)

print(head(conv$results))

tensor <- list(
  features_sig = sig_geom$features,
  features_int = reg_geom$features,
  meta         = sig_geom$meta
)

embedding <- compute_pca_embedding_from_tensor(tensor, n_components = 3)

# qual circuitry você quer?
cid <- "PAAD-3903 / PAAD-7083"
index <- which(rownames(sig_geom$features) == cid)

if (length(index) != 1L) {
  stop("Circuitries_id não encontrado ou não-único: ", cid,
       " | matches = ", length(index))
}

poly <- build_circuitry_polytope(tensor, embedding, index = index)

plot_signature_polytope_3d(poly)
plot_regulator_polytope_3d(poly)
plot_circuitry_3d(poly)

cat("OK: GeoCircuitry rodou em 18D (paper-aligned) com colnames alinhadas.\n")
