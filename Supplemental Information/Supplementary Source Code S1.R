###############################################################################
# Supplementary Source Code S1
# SigPolytope — Reviewer-1 Validation Additions (Deterministic; No Simulation)
#
# PURPOSE:
#   Uses ONLY columns already present in Dataset_S1 to:
#   (i)   reconstruct 18D latent matrices for sig and int
#   (ii)  run PCA and recompute barycenter-distance (dbary) in PC1–PC3
#   (iii) run a permutation-based null model for dbary (shuffle int within CTAB)
#   (iv)  compute baseline (non-geometric) metrics and a correlation heatmap
#   (v)   perform robustness variants of the embedding/scaling choices
#   (vi)  write audit artifacts + Supplementary Table S1
#   (vii) bundle ALL outputs into ONE Excel workbook (Dataset S2) with README sheet
#
# INPUTS (must exist in working directory or provide full paths):
#   - Dataset_S1.tsv
#   - Regulatory_circuitries.tsv  (optional; imported but not used)
###############################################################################

suppressPackageStartupMessages({
  library(rio)
})

options(stringsAsFactors = FALSE)

# ---- Output directory (ALWAYS ensure it exists, even if running block-by-block)
OUTDIR <- "Reviewer_1_validation_outputs"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("R version:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n")

# ---- Session info (audit)
writeLines(
  c(
    paste("R version:", R.version.string),
    paste("Platform:", R.version$platform),
    paste("Running:", Sys.time())
  ),
  con = file.path(OUTDIR, "session_info.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTDIR, "sessionInfo_full.txt")
)

# ---- Inputs
Dataset_S1 <- rio::import("Dataset_S1.tsv")
Dataset_S1 <- as.data.frame(Dataset_S1)

# Optional cross-check dataset; not used below
Regulatory_circuitries <- rio::import("Regulatory_circuitries.tsv")

###############################################################################
# 0) Hard preflight: required columns
###############################################################################

required_cols <- c(
  # identifiers / stratification
  "Circuitries_id", "CTAB",
  
  # geometric outputs already in Dataset S1 (used for cross-check)
  "barycenter_distance", "distance_implication",
  "sig_hull_vol", "int_hull_vol", "vol_ratio", "mean_volume",
  
  # correlation axes (rho + adjusted p)
  "Correlation_rho_sig", "Correlation_p.adj_sig",
  "Correlation_rho_int", "Correlation_p.adj_int",
  
  # tumor-normal axes (direction + adjusted p)
  "Tumor_vs_normal_sig", "Tumor_vs_normal_p.adj_sig",
  "Tumor_vs_normal_int", "Tumor_vs_normal_p.adj_int",
  
  # Cox directions + p-values (4 endpoints)
  "Cox_OS_type_sig",  "Cox_OS_p.value_sig",
  "Cox_OS_type_int",  "Cox_OS_p.value_int",
  "Cox_DSS_type_sig", "Cox_DSS_p.value_sig",
  "Cox_DSS_type_int", "Cox_DSS_p.value_int",
  "Cox_DFI_type_sig", "Cox_DFI_p.value_sig",
  "Cox_DFI_type_int", "Cox_DFI_p.value_int",
  "Cox_PFI_type_sig", "Cox_PFI_p.value_sig",
  "Cox_PFI_type_int", "Cox_PFI_p.value_int",
  
  # log-rank chi-square per endpoint (per side)
  "OS_log_rank_chisq_sig",  "OS_log_rank_chisq_int",
  "DSS_log_rank_chisq_sig", "DSS_log_rank_chisq_int",
  "DFI_log_rank_chisq_sig", "DFI_log_rank_chisq_int",
  "PFI_log_rank_chisq_sig", "PFI_log_rank_chisq_int",
  
  # microenvironment score
  "Microenvironment_score_sig", "Microenvironment_score_int",
  
  # immune class per side
  "Immune_classification_sig", "Immune_classification_int"
)

missing_cols <- setdiff(required_cols, names(Dataset_S1))
if (length(missing_cols) > 0) {
  stop(
    "Dataset_S1 is missing required columns:\n- ",
    paste(missing_cols, collapse = "\n- ")
  )
}

###############################################################################
# 1) Deterministic encoding helpers
###############################################################################

as_num <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

p_to_strength <- function(p) {
  p <- as_num(p)
  out <- rep(0, length(p))
  ok <- is.finite(p) & !is.na(p) & (p > 0)
  out[ok] <- -log10(p[ok])
  out[!is.finite(out)] <- 0
  out
}

map_cox_dir <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  allowed <- c("Risky", "Protective", "NS", "No data")
  bad <- setdiff(unique(x), allowed)
  if (length(bad) > 0) {
    stop(
      "Unrecognized Cox_*_type label(s) in Dataset_S1: ",
      paste(bad, collapse = ", "),
      "\nAllowed labels: ",
      paste(allowed, collapse = ", ")
    )
  }
  out <- numeric(length(x))
  out[x == "Risky"]      <-  1
  out[x == "Protective"] <- -1
  out[x == "NS"]         <-  0
  out[x == "No data"]    <-  0
  out
}

map_tn_dir <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  allowed <- c("Overexpression", "Underexpression", "Unchanged", "No data")
  bad <- setdiff(unique(x), allowed)
  if (length(bad) > 0) {
    stop(
      "Unrecognized Tumor_vs_normal label(s) in Dataset_S1: ",
      paste(bad, collapse = ", "),
      "\nAllowed labels: ",
      paste(allowed, collapse = ", ")
    )
  }
  out <- numeric(length(x))
  out[x == "Overexpression"]  <-  1
  out[x == "Underexpression"] <- -1
  out[x == "Unchanged"]       <-  0
  out[x == "No data"]         <-  0
  out
}

map_immune_dir <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  allowed <- c("Hot", "Cold", "Variable")
  bad <- setdiff(unique(x), allowed)
  if (length(bad) > 0) {
    stop(
      "Unrecognized Immune_classification label(s) in Dataset_S1: ",
      paste(bad, collapse = ", "),
      "\nAllowed labels: ",
      paste(allowed, collapse = ", ")
    )
  }
  out <- numeric(length(x))
  out[x == "Hot"]      <-  1
  out[x == "Cold"]     <- -1
  out[x == "Variable"] <-  0
  out
}

neutralize_survival_axes <- function(dir, strength, chisq) {
  dir <- as_num(dir)
  strength <- as_num(strength)
  chisq <- as_num(chisq)
  
  dir[!is.finite(dir) | is.na(dir)] <- 0
  strength[!is.finite(strength) | is.na(strength)] <- 0
  chisq[!is.finite(chisq) | is.na(chisq)] <- 0
  
  idx0 <- (dir == 0)
  strength[idx0] <- 0
  chisq[idx0] <- 0
  
  list(dir = dir, strength = strength, chisq = chisq)
}

###############################################################################
# 2) Reconstruct 18D latent vectors for sig and int (ONLY from Dataset_S1)
###############################################################################

build_latent18 <- function(df, side = c("sig", "int")) {
  side <- match.arg(side)
  col <- function(stem) paste0(stem, "_", side)
  
  # Correlation axes
  rho <- as_num(df[[col("Correlation_rho")]])
  rho[!is.finite(rho) | is.na(rho)] <- 0
  rho_strength <- p_to_strength(df[[col("Correlation_p.adj")]])
  
  # Tumor-normal axes
  tn_dir <- map_tn_dir(df[[col("Tumor_vs_normal")]])
  tn_strength <- p_to_strength(df[[col("Tumor_vs_normal_p.adj")]])
  
  # Survival axes (4 endpoints × {dir, strength, chisq})
  build_endpoint <- function(ep) {
    dir <- map_cox_dir(df[[col(paste0("Cox_", ep, "_type"))]])
    strength <- p_to_strength(df[[col(paste0("Cox_", ep, "_p.value"))]])
    chisq <- as_num(df[[paste0(ep, "_log_rank_chisq_", side)]])
    neutralize_survival_axes(dir, strength, chisq)
  }
  
  OS  <- build_endpoint("OS")
  DSS <- build_endpoint("DSS")
  DFI <- build_endpoint("DFI")
  PFI <- build_endpoint("PFI")
  
  # Context axes
  tme <- as_num(df[[col("Microenvironment_score")]])
  tme[!is.finite(tme) | is.na(tme)] <- 0
  immune_dir <- map_immune_dir(df[[col("Immune_classification")]])
  
  X <- cbind(
    rho_dir      = rho,
    rho_strength = rho_strength,
    TN_dir       = tn_dir,
    TN_strength  = tn_strength,
    
    OS_dir       = OS$dir,
    OS_strength  = OS$strength,
    OS_chisq     = OS$chisq,
    
    DSS_dir      = DSS$dir,
    DSS_strength = DSS$strength,
    DSS_chisq    = DSS$chisq,
    
    DFI_dir      = DFI$dir,
    DFI_strength = DFI$strength,
    DFI_chisq    = DFI$chisq,
    
    PFI_dir      = PFI$dir,
    PFI_strength = PFI$strength,
    PFI_chisq    = PFI$chisq,
    
    TME_score    = tme,
    Immune_dir   = immune_dir
  )
  
  X[!is.finite(X) | is.na(X)] <- 0
  X
}

X_sig18 <- build_latent18(Dataset_S1, "sig")
X_int18 <- build_latent18(Dataset_S1, "int")

stopifnot(nrow(X_sig18) == nrow(Dataset_S1))
stopifnot(nrow(X_int18) == nrow(Dataset_S1))
stopifnot(ncol(X_sig18) == 18)
stopifnot(ncol(X_int18) == 18)

###############################################################################
# 3) PCA embedding + dbary recomputation (PC1–PC3 distance)
###############################################################################

compute_dbary_from_18D <- function(X_sig, X_int, center = TRUE, scale. = TRUE) {
  X_all <- rbind(X_sig, X_int)
  pca <- prcomp(X_all, center = center, scale. = scale.)
  scores <- pca$x[, 1:3, drop = FALSE]
  n <- nrow(X_sig)
  sig_pc <- scores[seq_len(n), , drop = FALSE]
  int_pc <- scores[n + seq_len(n), , drop = FALSE]
  dbary <- sqrt(rowSums((sig_pc - int_pc)^2))
  list(pca = pca, sig_pc = sig_pc, int_pc = int_pc, dbary = dbary)
}

geom_main <- compute_dbary_from_18D(X_sig18, X_int18, center = TRUE, scale. = TRUE)
dbary_recomputed <- geom_main$dbary

dbary_original <- as_num(Dataset_S1$barycenter_distance)
ok_idx <- is.finite(dbary_original) & !is.na(dbary_original)
if (sum(ok_idx) < 10) stop("Not enough finite barycenter_distance values for cross-check.")

rho_check <- suppressWarnings(cor(dbary_original[ok_idx], dbary_recomputed[ok_idx], method = "spearman"))
cat(sprintf("Spearman correlation (Dataset_S1 barycenter_distance vs recomputed dbary_pca): %.4f\n", rho_check))

###############################################################################
# 4) Baseline metrics + correlation heatmap
###############################################################################

euclid18 <- sqrt(rowSums((X_sig18 - X_int18)^2))

cosine_dist <- function(A, B) {
  dot <- rowSums(A * B)
  na <- sqrt(rowSums(A * A))
  nb <- sqrt(rowSums(B * B))
  denom <- na * nb
  out <- rep(0, length(dot))
  ok <- (denom > 0) & is.finite(denom)
  out[ok] <- 1 - (dot[ok] / denom[ok])
  out[!is.finite(out) | is.na(out)] <- 0
  out
}
cos18 <- cosine_dist(X_sig18, X_int18)

immune_mismatch <- as.integer(
  trimws(tolower(as.character(Dataset_S1$Immune_classification_sig))) !=
    trimws(tolower(as.character(Dataset_S1$Immune_classification_int)))
)

tn_mismatch <- as.integer(
  trimws(tolower(as.character(Dataset_S1$Tumor_vs_normal_sig))) !=
    trimws(tolower(as.character(Dataset_S1$Tumor_vs_normal_int)))
)

cox_dir_sig <- cbind(
  OS  = map_cox_dir(Dataset_S1$Cox_OS_type_sig),
  DSS = map_cox_dir(Dataset_S1$Cox_DSS_type_sig),
  DFI = map_cox_dir(Dataset_S1$Cox_DFI_type_sig),
  PFI = map_cox_dir(Dataset_S1$Cox_PFI_type_sig)
)
cox_dir_int <- cbind(
  OS  = map_cox_dir(Dataset_S1$Cox_OS_type_int),
  DSS = map_cox_dir(Dataset_S1$Cox_DSS_type_int),
  DFI = map_cox_dir(Dataset_S1$Cox_DFI_type_int),
  PFI = map_cox_dir(Dataset_S1$Cox_PFI_type_int)
)
surv_dir_mismatch_count <- rowSums(sign(cox_dir_sig) != sign(cox_dir_int))

tme_abs_diff <- abs(as_num(Dataset_S1$Microenvironment_score_sig) - as_num(Dataset_S1$Microenvironment_score_int))
tme_abs_diff[!is.finite(tme_abs_diff) | is.na(tme_abs_diff)] <- 0

mean_vol <- as_num(Dataset_S1$mean_volume)
vol_ratio <- as_num(Dataset_S1$vol_ratio)

metrics_df <- data.frame(
  dbary_pca = dbary_recomputed,
  euclid18  = euclid18,
  cos18     = cos18,
  surv_dir_mismatch_count = surv_dir_mismatch_count,
  immune_mismatch = immune_mismatch,
  tn_mismatch = tn_mismatch,
  tme_abs_diff = tme_abs_diff,
  mean_volume = mean_vol,
  vol_ratio = vol_ratio,
  stringsAsFactors = FALSE
)

cor_mat <- suppressWarnings(cor(metrics_df, method = "spearman", use = "pairwise.complete.obs"))

rho_dbary_euclid <- cor_mat["dbary_pca", "euclid18"]
rho_dbary_cos    <- cor_mat["dbary_pca", "cos18"]
rho_dbary_surv   <- cor_mat["dbary_pca", "surv_dir_mismatch_count"]

cat(sprintf("Spearman rho(dbary, Euclid18)  = %.3f\n", rho_dbary_euclid))
cat(sprintf("Spearman rho(dbary, Cosine18)  = %.3f\n", rho_dbary_cos))
cat(sprintf("Spearman rho(dbary, SurvMismatch)= %.3f\n", rho_dbary_surv))

plot_cor_heatmap <- function(cor_mat,
                             main = "Spearman correlations: geometric vs baseline metrics",
                             cex_row = 0.95,
                             cex_col = 0.95,
                             cex_main = 1.15) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(12, 12, 6, 14), xpd = NA, cex.main = cex_main)
  heatmap(
    cor_mat,
    Rowv = NA, Colv = NA, scale = "none",
    main = main,
    margins = c(12, 12),
    cexRow = cex_row, cexCol = cex_col
  )
}

tiff(
  filename = file.path(OUTDIR, "metrics_spearman_heatmap.tiff"),
  width = 7000, height = 5000, res = 600,
  compression = "lzw"
)
plot_cor_heatmap(cor_mat)
dev.off()

pdf(
  file = file.path(OUTDIR, "metrics_spearman_heatmap.pdf"),
  width = 9.5, height = 7.0,
  useDingbats = FALSE
)
plot_cor_heatmap(cor_mat)
dev.off()

###############################################################################
# 5) Permutation null model for dbary (shuffle int within CTAB)
###############################################################################

classify_dbary_regime <- function(db) {
  cut(
    db,
    breaks = c(-Inf, 0.5, 1.5, 2.5, Inf),
    labels = c("High_concordance", "Moderate_discordance", "Strong_discordance", "Extreme_discordance"),
    right = FALSE
  )
}

observed_regime <- classify_dbary_regime(dbary_recomputed)
observed_prop <- prop.table(table(observed_regime))

permute_dbary_once <- function(X_sig, X_int, ctab, center = TRUE, scale. = TRUE) {
  idx <- seq_along(ctab)
  perm_idx <- idx
  for (g in unique(ctab)) {
    ii <- which(ctab == g)
    if (length(ii) > 1) perm_idx[ii] <- sample(ii, size = length(ii), replace = FALSE)
  }
  X_int_perm <- X_int[perm_idx, , drop = FALSE]
  geom <- compute_dbary_from_18D(X_sig, X_int_perm, center = center, scale. = scale.)
  geom$dbary
}

set.seed(1)

B <- 500
ctab <- as.character(Dataset_S1$CTAB)

null_mean_dbary   <- numeric(B)
null_prop_extreme <- numeric(B)
null_prop_high    <- numeric(B)

for (b in seq_len(B)) {
  db_null <- permute_dbary_once(X_sig18, X_int18, ctab, center = TRUE, scale. = TRUE)
  null_mean_dbary[b] <- mean(db_null, na.rm = TRUE)
  
  reg_null <- classify_dbary_regime(db_null)
  tab <- prop.table(table(reg_null))
  
  null_prop_extreme[b] <- if ("Extreme_discordance" %in% names(tab)) tab[["Extreme_discordance"]] else 0
  null_prop_high[b]    <- if ("High_concordance" %in% names(tab)) tab[["High_concordance"]] else 0
  
  if (b %% 50 == 0) cat("Permutation replicate:", b, "of", B, "\n")
}

obs_mean_dbary <- mean(dbary_recomputed, na.rm = TRUE)

get_prop <- function(prop_table, name) {
  if (!name %in% names(prop_table)) return(0)
  as.numeric(prop_table[[name]])
}

obs_extreme <- get_prop(observed_prop, "Extreme_discordance")
obs_high    <- get_prop(observed_prop, "High_concordance")

p_mean_dbary <- (sum(null_mean_dbary   >= obs_mean_dbary) + 1) / (B + 1)
p_extreme    <- (sum(null_prop_extreme >= obs_extreme) + 1) / (B + 1)
p_high       <- (sum(null_prop_high    >= obs_high) + 1) / (B + 1)

cat(sprintf("Observed mean dbary: %.4f | Empirical p (null >= obs): %.4g\n", obs_mean_dbary, p_mean_dbary))
cat(sprintf("Observed extreme proportion: %.4f | Empirical p (null >= obs): %.4g\n", obs_extreme, p_extreme))
cat(sprintf("Observed high-concordance proportion: %.4f | Empirical p (null >= obs): %.4g\n", obs_high, p_high))

tiff(
  filename = file.path(OUTDIR, "dbary_permutation_diagnostics.tiff"),
  width = 2400, height = 1200, res = 600,
  compression = "lzw"
)
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))

hist(
  null_mean_dbary, breaks = 40,
  main = "Null distribution of mean(dbary)",
  xlab = "mean(dbary) under CTAB-shuffled int",
  border = "white"
)
abline(v = obs_mean_dbary, lwd = 2)

hist(
  null_prop_extreme, breaks = 40,
  main = "Null distribution of P(Extreme)",
  xlab = "Extreme_discordance fraction under null",
  border = "white"
)
abline(v = obs_extreme, lwd = 2)

par(op)
dev.off()

###############################################################################
# 6) Robustness variants (embedding/scaling)
###############################################################################

compute_robustness <- function(X_sig, X_int, variant_name,
                               center = TRUE, scale. = TRUE,
                               drop_cols = character(0),
                               dbary_ref = dbary_recomputed) {
  if (length(drop_cols) > 0) {
    stopifnot(all(drop_cols %in% colnames(X_sig)))
    keep <- setdiff(colnames(X_sig), drop_cols)
    X_sig2 <- X_sig[, keep, drop = FALSE]
    X_int2 <- X_int[, keep, drop = FALSE]
  } else {
    X_sig2 <- X_sig
    X_int2 <- X_int
  }
  
  geom <- compute_dbary_from_18D(X_sig2, X_int2, center = center, scale. = scale.)
  db <- geom$dbary
  
  ok <- is.finite(db) & is.finite(dbary_ref)
  rho <- suppressWarnings(cor(db[ok], dbary_ref[ok], method = "spearman"))
  
  reg_ref <- classify_dbary_regime(dbary_ref)
  reg_new <- classify_dbary_regime(db)
  agree <- mean(as.character(reg_ref) == as.character(reg_new), na.rm = TRUE)
  
  list(variant = variant_name, spearman_rho = rho, regime_agreement = agree, dbary = db)
}

rob1 <- compute_robustness(
  X_sig18, X_int18,
  variant_name = "PCA_no_scaling",
  center = TRUE, scale. = FALSE
)

survival_cols <- c(
  "OS_dir","OS_strength","OS_chisq",
  "DSS_dir","DSS_strength","DSS_chisq",
  "DFI_dir","DFI_strength","DFI_chisq",
  "PFI_dir","PFI_strength","PFI_chisq"
)

rob2 <- compute_robustness(
  X_sig18, X_int18,
  variant_name = "Drop_survival_axes",
  center = TRUE, scale. = TRUE,
  drop_cols = survival_cols
)

rob3 <- compute_robustness(
  X_sig18, X_int18,
  variant_name = "Drop_immune_axis",
  center = TRUE, scale. = TRUE,
  drop_cols = c("Immune_dir")
)

robustness_summary <- data.frame(
  variant = c(rob1$variant, rob2$variant, rob3$variant),
  spearman_rho = c(rob1$spearman_rho, rob2$spearman_rho, rob3$spearman_rho),
  regime_agreement = c(rob1$regime_agreement, rob2$regime_agreement, rob3$regime_agreement),
  stringsAsFactors = FALSE
)
print(robustness_summary)

###############################################################################
# 7) Write audit artifacts (CSV) + null summaries
###############################################################################

latent_sig_out <- data.frame(Circuitries_id = Dataset_S1$Circuitries_id, X_sig18, check.names = FALSE)
latent_int_out <- data.frame(Circuitries_id = Dataset_S1$Circuitries_id, X_int18, check.names = FALSE)

write.csv(latent_sig_out, file = file.path(OUTDIR, "latent18_sig_reconstructed.csv"), row.names = FALSE)
write.csv(latent_int_out, file = file.path(OUTDIR, "latent18_int_reconstructed.csv"), row.names = FALSE)

write.csv(metrics_df, file = file.path(OUTDIR, "baseline_and_geometric_metrics.csv"), row.names = FALSE)

write.csv(as.data.frame(cor_mat),
          file = file.path(OUTDIR, "metrics_spearman_correlation_matrix.csv"),
          row.names = TRUE)

write.csv(robustness_summary,
          file = file.path(OUTDIR, "embedding_robustness_summary.csv"),
          row.names = FALSE)

null_summary <- data.frame(
  B = B,
  observed_mean_dbary = obs_mean_dbary,
  p_mean_dbary = p_mean_dbary,
  observed_extreme_prop = obs_extreme,
  p_extreme_prop = p_extreme,
  observed_high_concordance_prop = obs_high,
  p_high_concordance_prop = p_high,
  null_mean_dbary_mean = mean(null_mean_dbary),
  null_mean_dbary_sd = sd(null_mean_dbary),
  null_extreme_prop_mean = mean(null_prop_extreme),
  null_extreme_prop_sd = sd(null_prop_extreme),
  null_high_prop_mean = mean(null_prop_high),
  null_high_prop_sd = sd(null_prop_high),
  stringsAsFactors = FALSE
)

write.csv(null_summary,
          file = file.path(OUTDIR, "dbary_permutation_null_summary.csv"),
          row.names = FALSE)

write.csv(
  data.frame(
    null_mean_dbary   = null_mean_dbary,
    null_prop_extreme = null_prop_extreme,
    null_prop_high    = null_prop_high
  ),
  file = file.path(OUTDIR, "dbary_permutation_null_raw.csv"),
  row.names = FALSE
)

###############################################################################
# 8) Supplementary Table S1 (summary table) + Dataset S2 (single Excel workbook)
###############################################################################

# ---- Interpretation helper for Table S1
interpret_rho <- function(r) {
  ar <- abs(as.numeric(r))
  if (!is.finite(ar) || is.na(ar)) return(NA_character_)
  if (ar < 0.10) return("negligible")
  if (ar < 0.30) return("weak")
  if (ar < 0.50) return("moderate")
  if (ar < 0.70) return("strong")
  "very strong"
}

# ---- Table S1 core correlations
stopifnot("dbary_pca" %in% rownames(cor_mat))

required_metrics <- c("euclid18","cos18","surv_dir_mismatch_count","immune_mismatch","tn_mismatch","tme_abs_diff")
missing_m <- setdiff(required_metrics, colnames(cor_mat))
if (length(missing_m) > 0) stop("cor_mat missing required metric(s): ", paste(missing_m, collapse = ", "))

S1_core <- data.frame(
  Comparison = c(
    "dbary (PC1–PC3) vs Euclidean distance (latent 18D)",
    "dbary (PC1–PC3) vs cosine distance (latent 18D)",
    "dbary (PC1–PC3) vs survival-direction mismatch count (OS/DSS/DFI/PFI)",
    "dbary (PC1–PC3) vs immune classification mismatch",
    "dbary (PC1–PC3) vs tumor–normal direction mismatch",
    "dbary (PC1–PC3) vs microenvironment score absolute difference"
  ),
  Metric_A = rep("dbary_pca", 6),
  Metric_B = c("euclid18","cos18","surv_dir_mismatch_count","immune_mismatch","tn_mismatch","tme_abs_diff"),
  Spearman_rho = c(
    cor_mat["dbary_pca","euclid18"],
    cor_mat["dbary_pca","cos18"],
    cor_mat["dbary_pca","surv_dir_mismatch_count"],
    cor_mat["dbary_pca","immune_mismatch"],
    cor_mat["dbary_pca","tn_mismatch"],
    cor_mat["dbary_pca","tme_abs_diff"]
  ),
  N = nrow(metrics_df),
  stringsAsFactors = FALSE
)
S1_core$Interpretation <- vapply(S1_core$Spearman_rho, interpret_rho, character(1))

S1_repro <- data.frame(
  Comparison = "Reproducibility cross-check: Dataset_S1 barycenter_distance vs recomputed dbary_pca",
  Metric_A = "Dataset_S1$barycenter_distance",
  Metric_B = "dbary_pca (recomputed)",
  Spearman_rho = as.numeric(rho_check),
  N = sum(ok_idx),
  Interpretation = interpret_rho(rho_check),
  stringsAsFactors = FALSE
)

S1_robust <- data.frame(
  Comparison = paste0("Embedding robustness variant: ", robustness_summary$variant),
  Metric_A = "dbary_pca (main)",
  Metric_B = "dbary_pca (variant)",
  Spearman_rho = as.numeric(robustness_summary$spearman_rho),
  N = nrow(metrics_df),
  Interpretation = vapply(as.numeric(robustness_summary$spearman_rho), interpret_rho, character(1)),
  stringsAsFactors = FALSE
)
S1_robust$Regime_agreement <- as.numeric(robustness_summary$regime_agreement)

S1_perm <- data.frame(
  Comparison = "Permutation null (CTAB-shuffled int): mean(dbary) and regime proportions",
  Metric_A = "dbary_pca",
  Metric_B = "CTAB-shuffled int null",
  Spearman_rho = NA_real_,
  N = nrow(metrics_df),
  Interpretation = NA_character_,
  Regime_agreement = NA_real_,
  Observed_mean_dbary = as.numeric(null_summary$observed_mean_dbary),
  Null_mean_dbary_mean = as.numeric(null_summary$null_mean_dbary_mean),
  P_value_mean_dbary = as.numeric(null_summary$p_mean_dbary),
  Observed_extreme_prop = as.numeric(null_summary$observed_extreme_prop),
  P_value_extreme_prop = as.numeric(null_summary$p_extreme_prop),
  B = as.numeric(null_summary$B),
  stringsAsFactors = FALSE
)

# ---- Schema-aligned stacking
align_cols <- function(df, all_cols) {
  miss <- setdiff(all_cols, names(df))
  if (length(miss) > 0) for (m in miss) df[[m]] <- NA
  df[, all_cols, drop = FALSE]
}
all_cols <- unique(c(names(S1_core), names(S1_repro), names(S1_robust), names(S1_perm)))

S1 <- rbind(
  align_cols(S1_core,  all_cols),
  align_cols(S1_repro, all_cols),
  align_cols(S1_robust, all_cols),
  align_cols(S1_perm,  all_cols)
)

front <- c("Comparison","Metric_A","Metric_B","Spearman_rho","N","Interpretation","Regime_agreement")
front <- intersect(front, names(S1))
S1 <- S1[, c(front, setdiff(names(S1), front)), drop = FALSE]

# ------------------------------------------------------------------
# Build Excel workbook (README FIRST)
# ------------------------------------------------------------------

wb <- openxlsx::createWorkbook()

# 1) README sheet (FIRST TAB)
openxlsx::addWorksheet(wb, "README")

sheet_map <- data.frame(
  Sheet = c(
    "SuppTable_S1",
    "baseline_geom_metrics",
    "spearman_corr_matrix",
    "embedding_robustness",
    "perm_null_summary",
    "perm_null_raw",
    "latent18_sig",
    "latent18_int"
  ),
  Description = c(
    "Programmatically generated summary table reporting benchmarking results comparing the geometric discordance metric (dbary) with baseline discordance measures, including reproducibility checks, robustness analyses, and permutation-based null statistics.",
    "Dataset containing all computed discordance metrics for each regulatory circuitry, including dbary, Euclidean and cosine distances in the reconstructed 18D latent space, survival-direction mismatches across OS/DSS/DFI/PFI, immune mismatches, tumor–normal mismatches, and microenvironment score differences.",
    "Spearman correlation matrix comparing the geometric discordance metric with all baseline discordance metrics.",
    "Embedding robustness analysis showing stability of dbary under alternative reconstruction choices such as PCA without scaling or removal of specific biological axes.",
    "Summary statistics from the permutation-based null model used to test whether observed geometric discordance exceeds expectations under random pairing.",
    "Raw outputs from the permutation null simulations (one row per permutation replicate).",
    "Reconstructed 18-dimensional latent vectors representing the phenotypic embedding of molecular signatures.",
    "Reconstructed 18-dimensional latent vectors representing the phenotypic embedding of regulatory interactions."
  ),
  stringsAsFactors = FALSE
)

openxlsx::writeData(wb, "README", sheet_map)

# ------------------------------------------------------------------
# Variable dictionary sheet (programmatic, human-readable)
# ------------------------------------------------------------------

# Helper: dictionary of variable descriptions
var_dictionary <- c(
  # ---- Supplementary Table S1
  Comparison = "Human-readable description of the comparison or validation test summarized in the row.",
  Metric_A = "Primary metric being evaluated in the comparison.",
  Metric_B = "Reference or alternative metric used for comparison.",
  Spearman_rho = "Spearman rank correlation coefficient quantifying monotonic association between two metrics.",
  N = "Number of circuitries included in the calculation.",
  Interpretation = "Qualitative interpretation of the correlation magnitude.",
  Regime_agreement = "Proportion of circuitries assigned to the same discordance regime under two alternative embeddings or reconstruction variants.",
  Observed_mean_dbary = "Observed mean barycenter-distance metric (dbary) across all circuitries.",
  Null_mean_dbary_mean = "Mean of the null distribution of the average dbary obtained from CTAB-stratified permutations.",
  P_value_mean_dbary = "Empirical p-value for the observed mean dbary relative to the permutation null distribution.",
  Observed_extreme_prop = "Observed proportion of circuitries classified in the extreme-discordance regime.",
  P_value_extreme_prop = "Empirical p-value for the observed extreme-discordance proportion relative to the permutation null distribution.",
  B = "Number of permutation replicates used in the null model.",
  
  # ---- baseline_geom_metrics
  dbary_pca = "Geometric discordance metric defined as the Euclidean distance between signature and interaction barycenters in the shared PC1–PC3 embedding.",
  euclid18 = "Euclidean distance between reconstructed signature and interaction vectors in the 18-dimensional latent space.",
  cos18 = "Cosine distance between reconstructed signature and interaction vectors in the 18-dimensional latent space.",
  surv_dir_mismatch_count = "Count of survival endpoints (OS, DSS, DFI, PFI) for which signature and interaction components show discordant survival-direction encoding.",
  immune_mismatch = "Indicator of mismatch between immune classifications assigned to signature and interaction components.",
  tn_mismatch = "Indicator of mismatch between tumor–normal direction labels assigned to signature and interaction components.",
  tme_abs_diff = "Absolute difference between microenvironment scores of signature and interaction components.",
  mean_volume = "Mean convex-hull volume across signature and interaction polytopes for a given circuitry.",
  vol_ratio = "Ratio between signature and interaction convex-hull volumes, used to quantify geometric asymmetry.",
  
  # ---- robustness
  variant = "Name of the robustness-analysis condition or embedding variant.",
  regime_agreement = "Fraction of circuitries assigned to the same discordance regime under the alternative variant and the reference embedding.",
  
  # ---- permutation summary/raw
  p_mean_dbary = "Empirical p-value for the observed mean dbary relative to the permutation null distribution.",
  observed_extreme_prop = "Observed proportion of circuitries classified as extreme discordance.",
  p_extreme_prop = "Empirical p-value for the observed extreme-discordance proportion.",
  observed_high_concordance_prop = "Observed proportion of circuitries classified as highly concordant.",
  p_high_concordance_prop = "Empirical p-value for the observed high-concordance proportion.",
  null_mean_dbary_mean = "Mean of the null distribution of average dbary values across permutation replicates.",
  null_mean_dbary_sd = "Standard deviation of the null distribution of average dbary values.",
  null_extreme_prop_mean = "Mean null proportion of circuitries classified as extreme discordance.",
  null_extreme_prop_sd = "Standard deviation of the null proportion of extreme discordance.",
  null_high_prop_mean = "Mean null proportion of circuitries classified as high concordance.",
  null_high_prop_sd = "Standard deviation of the null proportion of high concordance.",
  null_mean_dbary = "Mean dbary obtained in a single permutation replicate.",
  null_prop_extreme = "Proportion of circuitries classified as extreme discordance in a single permutation replicate.",
  null_prop_high = "Proportion of circuitries classified as high concordance in a single permutation replicate.",
  
  # ---- latent 18D vectors
  Circuitries_id = "Unique identifier of the regulatory circuitry.",
  rho_dir = "Signed correlation coefficient describing the direction and magnitude of association in the correlation axis.",
  rho_strength = "Strength of the correlation signal, encoded as -log10(adjusted p-value).",
  TN_dir = "Tumor–normal direction encoded as +1 (overexpression), -1 (underexpression), or 0 (unchanged/no data).",
  TN_strength = "Strength of the tumor–normal association, encoded as -log10(adjusted p-value).",
  OS_dir = "Direction of overall-survival association encoded as +1 (risky), -1 (protective), or 0 (non-significant/no data).",
  OS_strength = "Strength of the overall-survival association, encoded as -log10(p-value).",
  OS_chisq = "Log-rank chi-square statistic for overall survival.",
  DSS_dir = "Direction of disease-specific survival association encoded as +1 (risky), -1 (protective), or 0 (non-significant/no data).",
  DSS_strength = "Strength of the disease-specific survival association, encoded as -log10(p-value).",
  DSS_chisq = "Log-rank chi-square statistic for disease-specific survival.",
  DFI_dir = "Direction of disease-free interval association encoded as +1 (risky), -1 (protective), or 0 (non-significant/no data).",
  DFI_strength = "Strength of the disease-free interval association, encoded as -log10(p-value).",
  DFI_chisq = "Log-rank chi-square statistic for disease-free interval.",
  PFI_dir = "Direction of progression-free interval association encoded as +1 (risky), -1 (protective), or 0 (non-significant/no data).",
  PFI_strength = "Strength of the progression-free interval association, encoded as -log10(p-value).",
  PFI_chisq = "Log-rank chi-square statistic for progression-free interval.",
  TME_score = "Microenvironment score associated with the circuitry component.",
  Immune_dir = "Immune-classification axis encoded as +1 (hot), -1 (cold), or 0 (variable)."
)

# Objects that generate workbook sheets
sheet_objects <- list(
  SuppTable_S1 = S1,
  baseline_geom_metrics = metrics_df,
  spearman_corr_matrix = as.data.frame(cor_mat),
  embedding_robustness = robustness_summary,
  perm_null_summary = null_summary,
  perm_null_raw = data.frame(
    null_mean_dbary   = null_mean_dbary,
    null_prop_extreme = null_prop_extreme,
    null_prop_high    = null_prop_high
  ),
  latent18_sig = latent_sig_out,
  latent18_int = latent_int_out
)

# Build variable dictionary programmatically
variable_dictionary <- do.call(
  rbind,
  lapply(names(sheet_objects), function(sh) {
    vars <- colnames(sheet_objects[[sh]])
    data.frame(
      Sheet = sh,
      Variable = vars,
      Description = ifelse(
        vars %in% names(var_dictionary),
        unname(var_dictionary[vars]),
        "Description not yet defined; please review this variable."
      ),
      stringsAsFactors = FALSE
    )
  })
)

# Optional: flag undocumented variables
undocumented_vars <- variable_dictionary$Variable[variable_dictionary$Description == "Description not yet defined; please review this variable."]
if (length(undocumented_vars) > 0) {
  cat("Warning: some variables do not yet have dictionary descriptions:\n")
  print(unique(undocumented_vars))
}

# Add dictionary sheet to workbook
openxlsx::addWorksheet(wb, "VARIABLE_DICTIONARY")
openxlsx::writeData(wb, "VARIABLE_DICTIONARY", variable_dictionary)

# 2) Supplementary Table S1
openxlsx::addWorksheet(wb, "SuppTable_S1")
openxlsx::writeData(wb, "SuppTable_S1", S1)

# 3) Baseline metrics
openxlsx::addWorksheet(wb, "baseline_geom_metrics")
openxlsx::writeData(wb, "baseline_geom_metrics", metrics_df)

# 4) Correlation matrix
openxlsx::addWorksheet(wb, "spearman_corr_matrix")
openxlsx::writeData(wb, "spearman_corr_matrix", as.data.frame(cor_mat), rowNames = TRUE)

# 5) Robustness results
openxlsx::addWorksheet(wb, "embedding_robustness")
openxlsx::writeData(wb, "embedding_robustness", robustness_summary)

# 6) Permutation null summary
openxlsx::addWorksheet(wb, "perm_null_summary")
openxlsx::writeData(wb, "perm_null_summary", null_summary)

# 7) Permutation raw outputs
openxlsx::addWorksheet(wb, "perm_null_raw")
openxlsx::writeData(
  wb, "perm_null_raw",
  data.frame(
    null_mean_dbary   = null_mean_dbary,
    null_prop_extreme = null_prop_extreme,
    null_prop_high    = null_prop_high
  )
)

# 8) Latent vectors (signature)
openxlsx::addWorksheet(wb, "latent18_sig")
openxlsx::writeData(wb, "latent18_sig", latent_sig_out)

# 9) Latent vectors (interaction)
openxlsx::addWorksheet(wb, "latent18_int")
openxlsx::writeData(wb, "latent18_int", latent_int_out)

# ------------------------------------------------------------------
# Save workbook
# ------------------------------------------------------------------

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

openxlsx::saveWorkbook(
  wb,
  file = file.path(OUTDIR, "SigPolytope_validation_outputs_ALL.xlsx"),
  overwrite = TRUE
)

openxlsx::saveWorkbook(
  wb,
  file = file.path(OUTDIR, "Dataset_S2.xlsx"),
  overwrite = TRUE
)

###############################################################################
# Extract key statistics for manuscript reporting
# Includes unbiased permutation p-values (+1 correction)
###############################################################################

# Load permutation summary
null_summary <- read.csv(
  file.path("Reviewer_1_validation_outputs", "dbary_permutation_null_summary.csv"),
  stringsAsFactors = FALSE
)

# Load raw permutation outputs
perm_raw <- read.csv(
  file.path("Reviewer_1_validation_outputs", "dbary_permutation_null_raw.csv"),
  stringsAsFactors = FALSE
)

###############################################################################
# Observed statistics
###############################################################################

obs_mean_dbary <- null_summary$observed_mean_dbary
obs_extreme    <- null_summary$observed_extreme_prop
obs_high       <- null_summary$observed_high_concordance_prop

###############################################################################
# Raw permutation distributions
###############################################################################

null_mean_dbary   <- perm_raw$null_mean_dbary
null_prop_extreme <- perm_raw$null_prop_extreme
null_prop_high    <- perm_raw$null_prop_high

B <- length(null_mean_dbary)

###############################################################################
# Unbiased permutation p-values (+1 correction)
###############################################################################

p_mean_dbary <- (sum(null_mean_dbary >= obs_mean_dbary) + 1) / (B + 1)
p_extreme    <- (sum(null_prop_extreme >= obs_extreme) + 1) / (B + 1)
p_high       <- (sum(null_prop_high >= obs_high) + 1) / (B + 1)

###############################################################################
# Helper function for manuscript-style p-value formatting
###############################################################################

format_p <- function(p, digits = 3) {
  if (!is.finite(p) || is.na(p)) return(NA_character_)
  if (p < 10^(-digits)) {
    return(paste0("< ", format(10^(-digits), scientific = FALSE)))
  } else {
    return(sprintf(paste0("%.", digits, "f"), p))
  }
}

###############################################################################
# Print manuscript-ready values
###############################################################################

cat("\n--- Permutation test summary (bias-corrected) ---\n")

cat(sprintf(
  "Observed mean dbary: %.4f\nEmpirical p-value: %s\n\n",
  obs_mean_dbary,
  format_p(p_mean_dbary)
))

cat(sprintf(
  "Observed extreme-discordance proportion: %.4f\nEmpirical p-value: %s\n\n",
  obs_extreme,
  format_p(p_extreme)
))

cat(sprintf(
  "Observed high-concordance proportion: %.4f\nEmpirical p-value: %s\n\n",
  obs_high,
  format_p(p_high)
))

cat("Done.\nOutputs written to: ", normalizePath(OUTDIR), "\n", sep = "")
###############################################################################

###############################################################################
###############################################################################
# Figure S5 — Permutation diagnostics (publication-ready, 3 panels)
#   Panel A: null distribution of mean(dbary)
#   Panel B: null distribution of extreme discordance proportion
#   Panel C: null distribution of high concordance proportion
#
# Layout:
#   A and B on the top row
#   C spanning the full bottom row
#
# Output:
#   - TIFF, 600 dpi, A4 landscape
#   - PDF companion file
###############################################################################

plot_perm_hist <- function(x,
                           obs,
                           main_title,
                           xlab,
                           fill_col,
                           panel_tag,
                           p_value,
                           breaks = 35,
                           obs_digits = 4,
                           p_digits = 3,
                           right_pad_frac = 0.16,
                           ann_x_inset = 0.04,
                           ann_y_inset = 0.08) {
  # Histogram without plotting first, so we can compute clean limits
  h <- hist(x, breaks = breaks, plot = FALSE)
  
  # Add extra right-side padding so the observed line/text do not feel cramped
  xr <- range(c(h$breaks, obs), finite = TRUE)
  xpad <- right_pad_frac * diff(xr)
  if (!is.finite(xpad) || xpad <= 0) xpad <- 0.1
  
  ymax <- max(h$counts, na.rm = TRUE)
  ylim_top <- ymax * 1.18
  
  hist(
    x,
    breaks = breaks,
    freq = TRUE,
    col = fill_col,
    border = "white",
    main = main_title,
    xlab = xlab,
    ylab = "Frequency",
    cex.main = 1.05,
    cex.lab = 0.95,
    cex.axis = 0.90,
    las = 1,
    xlim = c(xr[1], xr[2] + xpad),
    ylim = c(0, ylim_top)
  )
  
  abline(v = obs, lwd = 2.4, lty = 1)
  
  # Panel tag
  usr <- par("usr")
  text(
    x = usr[1] + 0.015 * diff(usr[1:2]),
    y = usr[4] - 0.04 * diff(usr[3:4]),
    labels = panel_tag,
    adj = c(0, 1),
    cex = 1.15,
    font = 2
  )
  
  # Observed statistic + empirical p-value annotation
  ann_text <- paste0(
    "Observed = ", formatC(obs, format = "f", digits = obs_digits),
    "\nEmpirical p = ", formatC(p_value, format = "f", digits = p_digits)
  )
  
  text(
    x = usr[2] - ann_x_inset * diff(usr[1:2]),
    y = usr[4] - ann_y_inset * diff(usr[3:4]),
    labels = ann_text,
    adj = c(1, 1),
    cex = 0.88
  )
}

make_figure_S5 <- function(tiff_file,
                           pdf_file,
                           null_mean_dbary,
                           null_prop_extreme,
                           null_prop_high,
                           obs_mean_dbary,
                           obs_extreme,
                           obs_high,
                           p_mean_dbary,
                           p_extreme,
                           p_high) {
  
  # A4 landscape in inches
  a4_width_in  <- 11.69
  a4_height_in <- 8.27
  
  # Shared drawing routine
  draw_panels <- function() {
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    # Top row: A, B
    # Bottom row: C spanning both columns
    layout(
      matrix(c(1, 2,
               3, 3), nrow = 2, byrow = TRUE),
      heights = c(1, 1.05),
      widths  = c(1, 1)
    )
    
    par(
      mar = c(4.6, 4.8, 3.2, 1.2),
      oma = c(0.5, 0.5, 0.5, 0.5)
    )
    
    plot_perm_hist(
      x = null_mean_dbary,
      obs = obs_mean_dbary,
      main_title = "Null distribution of mean(dbary)",
      xlab = "Mean dbary under CTAB-constrained permutation",
      fill_col = "#4E79A7",   # blue
      panel_tag = "A",
      p_value = p_mean_dbary,
      breaks = 35,
      obs_digits = 4,
      p_digits = 3,
      right_pad_frac = 0.16,
      ann_x_inset = 0.04,
      ann_y_inset = 0.08
    )
    
    plot_perm_hist(
      x = null_prop_extreme,
      obs = obs_extreme,
      main_title = "Null distribution of extreme discordance",
      xlab = "Proportion of extreme discordance under null",
      fill_col = "#E15759",   # red
      panel_tag = "B",
      p_value = p_extreme,
      breaks = 35,
      obs_digits = 4,
      p_digits = 3,
      right_pad_frac = 0.20,
      ann_x_inset = 0.05,
      ann_y_inset = 0.08
    )
    
    plot_perm_hist(
      x = null_prop_high,
      obs = obs_high,
      main_title = "Null distribution of high concordance",
      xlab = "Proportion of high concordance under null",
      fill_col = "#59A14F",   # green
      panel_tag = "C",
      p_value = p_high,
      breaks = 35,
      obs_digits = 4,
      p_digits = 3,
      right_pad_frac = 0.16,
      ann_x_inset = 0.04,
      ann_y_inset = 0.08
    )
  }
  
  # TIFF export: A4 landscape, 600 dpi
  tiff(
    filename = tiff_file,
    width = a4_width_in,
    height = a4_height_in,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  draw_panels()
  dev.off()
  
  # PDF export for visual inspection / journal workflow
  pdf(
    file = pdf_file,
    width = a4_width_in,
    height = a4_height_in,
    useDingbats = FALSE
  )
  draw_panels()
  dev.off()
}

make_figure_S5(
  tiff_file = file.path(OUTDIR, "Figure_S5_permutation_diagnostics_A4_600dpi.tiff"),
  pdf_file  = file.path(OUTDIR, "Figure_S5_permutation_diagnostics_A4.pdf"),
  null_mean_dbary   = null_mean_dbary,
  null_prop_extreme = null_prop_extreme,
  null_prop_high    = null_prop_high,
  obs_mean_dbary    = obs_mean_dbary,
  obs_extreme       = obs_extreme,
  obs_high          = obs_high,
  p_mean_dbary      = p_mean_dbary,
  p_extreme         = p_extreme,
  p_high            = p_high
)

###############################################################################
# Figure_S4. Spearman correlation heatmap comparing geometric and baseline
# discordance metrics.
###############################################################################

euclid18 <- sqrt(rowSums((X_sig18 - X_int18)^2))

cosine_dist <- function(A, B) {
  dot <- rowSums(A * B)
  na <- sqrt(rowSums(A * A))
  nb <- sqrt(rowSums(B * B))
  denom <- na * nb
  out <- rep(0, length(dot))
  ok <- (denom > 0) & is.finite(denom)
  out[ok] <- 1 - (dot[ok] / denom[ok])
  out[!is.finite(out) | is.na(out)] <- 0
  out
}
cos18 <- cosine_dist(X_sig18, X_int18)

immune_mismatch <- as.integer(
  trimws(tolower(as.character(Dataset_S1$Immune_classification_sig))) !=
    trimws(tolower(as.character(Dataset_S1$Immune_classification_int)))
)

tn_mismatch <- as.integer(
  trimws(tolower(as.character(Dataset_S1$Tumor_vs_normal_sig))) !=
    trimws(tolower(as.character(Dataset_S1$Tumor_vs_normal_int)))
)

cox_dir_sig <- cbind(
  OS  = map_cox_dir(Dataset_S1$Cox_OS_type_sig),
  DSS = map_cox_dir(Dataset_S1$Cox_DSS_type_sig),
  DFI = map_cox_dir(Dataset_S1$Cox_DFI_type_sig),
  PFI = map_cox_dir(Dataset_S1$Cox_PFI_type_sig)
)

cox_dir_int <- cbind(
  OS  = map_cox_dir(Dataset_S1$Cox_OS_type_int),
  DSS = map_cox_dir(Dataset_S1$Cox_DSS_type_int),
  DFI = map_cox_dir(Dataset_S1$Cox_DFI_type_int),
  PFI = map_cox_dir(Dataset_S1$Cox_PFI_type_int)
)

surv_dir_mismatch_count <- rowSums(sign(cox_dir_sig) != sign(cox_dir_int))

tme_abs_diff <- abs(
  as_num(Dataset_S1$Microenvironment_score_sig) -
    as_num(Dataset_S1$Microenvironment_score_int)
)
tme_abs_diff[!is.finite(tme_abs_diff) | is.na(tme_abs_diff)] <- 0

mean_vol  <- as_num(Dataset_S1$mean_volume)
vol_ratio <- as_num(Dataset_S1$vol_ratio)

metrics_df <- data.frame(
  dbary_pca = dbary_recomputed,
  euclid18  = euclid18,
  cos18     = cos18,
  surv_dir_mismatch_count = surv_dir_mismatch_count,
  immune_mismatch = immune_mismatch,
  tn_mismatch = tn_mismatch,
  tme_abs_diff = tme_abs_diff,
  mean_volume = mean_vol,
  vol_ratio = vol_ratio,
  stringsAsFactors = FALSE
)

cor_mat <- suppressWarnings(
  cor(metrics_df, method = "spearman", use = "pairwise.complete.obs")
)

rho_dbary_euclid <- cor_mat["dbary_pca", "euclid18"]
rho_dbary_cos    <- cor_mat["dbary_pca", "cos18"]
rho_dbary_surv   <- cor_mat["dbary_pca", "surv_dir_mismatch_count"]

cat(sprintf("Spearman rho(dbary, Euclid18)       = %.3f\n", rho_dbary_euclid))
cat(sprintf("Spearman rho(dbary, Cosine18)       = %.3f\n", rho_dbary_cos))
cat(sprintf("Spearman rho(dbary, SurvMismatch)   = %.3f\n", rho_dbary_surv))

# ------------------------------------------------------------------
# Publication labels
# ------------------------------------------------------------------
pretty_labels <- c(
  dbary_pca = "dbary (PC1–PC3)",
  euclid18 = "Euclidean distance (18D)",
  cos18 = "Cosine distance (18D)",
  surv_dir_mismatch_count = "Survival-direction mismatch count",
  immune_mismatch = "Immune classification mismatch",
  tn_mismatch = "Tumor–normal mismatch",
  tme_abs_diff = "Microenvironment absolute difference",
  mean_volume = "Mean volume",
  vol_ratio = "Volume ratio"
)

# ------------------------------------------------------------------
# Helper: apply pretty labels
# ------------------------------------------------------------------
apply_label_map_to_matrix <- function(mat, label_map = NULL) {
  mat <- as.matrix(mat)
  
  if (!is.null(label_map)) {
    rn <- rownames(mat)
    cn <- colnames(mat)
    
    rownames(mat) <- ifelse(
      rn %in% names(label_map),
      unname(label_map[rn]),
      rn
    )
    
    colnames(mat) <- ifelse(
      cn %in% names(label_map),
      unname(label_map[cn]),
      cn
    )
  }
  
  mat
}

# ------------------------------------------------------------------
# Helper: convert rendered text extent (inches) to margin lines
# mar is expressed in lines; conversion uses csi * mex
# ------------------------------------------------------------------
inches_to_mar_lines <- function(x_in) {
  x_in / (par("csi") * par("mex"))
}

# ------------------------------------------------------------------
# Helper: compute dynamic margins from actual rendered axis labels
# ------------------------------------------------------------------
compute_dynamic_heatmap_margins <- function(
    x_labels,
    y_labels,
    cex_axis = 0.96,
    bottom_pad_lines = 3.5,
    left_pad_lines = 3.5,
    min_bottom = 10.0,
    min_left = 12.0,
    max_bottom = 24.0,
    max_left = 28.0
) {
  x_labels <- as.character(x_labels)
  y_labels <- as.character(y_labels)
  
  x_labels[is.na(x_labels)] <- ""
  y_labels[is.na(y_labels)] <- ""
  
  # For las = 2 axis labels, the label extent that matters most is the
  # rendered string width, because the text is rotated perpendicular
  # to the axis direction.
  x_extent_in <- if (length(x_labels)) {
    max(strwidth(x_labels, units = "inches", cex = cex_axis), na.rm = TRUE)
  } else {
    0
  }
  
  y_extent_in <- if (length(y_labels)) {
    max(strwidth(y_labels, units = "inches", cex = cex_axis), na.rm = TRUE)
  } else {
    0
  }
  
  # Convert to margin lines and add padding for tick-label breathing room
  bottom_lines <- inches_to_mar_lines(x_extent_in) + bottom_pad_lines
  left_lines   <- inches_to_mar_lines(y_extent_in) + left_pad_lines
  
  bottom_lines <- max(min_bottom, min(bottom_lines, max_bottom))
  left_lines   <- max(min_left,   min(left_lines,   max_left))
  
  list(
    bottom = bottom_lines,
    left   = left_lines
  )
}

# ------------------------------------------------------------------
# Publication-ready correlation heatmap with automatic margin scaling
# and no internal plot title
# ------------------------------------------------------------------
plot_cor_heatmap_pub <- function(
    cor_mat,
    label_map = NULL,
    n_colors = 201,
    cex_axis = 0.96,
    cex_cell = 0.82,
    key_label = "Spearman rho",
    key_axis_cex = 0.90,
    layout_widths = c(7.4, 1.25),
    heatmap_right_mar = 2.0,
    heatmap_top_mar = 1.5,
    key_left_mar = 1.5,
    key_right_mar = 6.0,
    key_bottom_mar = NULL,
    debug_margins = FALSE
) {
  stopifnot(is.matrix(cor_mat) || is.data.frame(cor_mat))
  
  mat <- apply_label_map_to_matrix(cor_mat, label_map = label_map)
  n <- nrow(mat)
  
  # Symmetric palette around zero
  zlim   <- c(-1, 1)
  breaks <- seq(zlim[1], zlim[2], length.out = n_colors + 1)
  cols   <- grDevices::colorRampPalette(c("#3B4CC0", "#FFFFFF", "#B40426"))(n_colors)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  # No outer title, so no outer margin is needed
  par(oma = c(0, 0, 0, 0))
  
  # Compute dynamic margins after cex choices are known
  dyn_mar <- compute_dynamic_heatmap_margins(
    x_labels = colnames(mat),
    y_labels = rownames(mat),
    cex_axis = cex_axis,
    bottom_pad_lines = 3.5,
    left_pad_lines = 3.5,
    min_bottom = 10.0,
    min_left = 12.0,
    max_bottom = 24.0,
    max_left = 28.0
  )
  
  # Harmonize key bottom margin with heatmap bottom margin so both panels align
  if (is.null(key_bottom_mar)) {
    key_bottom_mar <- dyn_mar$bottom
  }
  
  if (isTRUE(debug_margins)) {
    message(sprintf(
      "Dynamic margins computed: bottom = %.2f lines; left = %.2f lines",
      dyn_mar$bottom, dyn_mar$left
    ))
  }
  
  # Two-panel layout: main heatmap + color key
  layout(
    mat = matrix(c(1, 2), nrow = 1),
    widths = layout_widths
  )
  
  # ----------------------------------------------------------------
  # Main heatmap panel
  # ----------------------------------------------------------------
  par(
    mar = c(dyn_mar$bottom, dyn_mar$left, heatmap_top_mar, heatmap_right_mar),
    xpd = NA
  )
  
  image(
    x = seq_len(n),
    y = seq_len(n),
    z = t(mat[n:1, , drop = FALSE]),
    col = cols,
    breaks = breaks,
    zlim = zlim,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  
  abline(
    h = seq(0.5, n + 0.5, by = 1),
    v = seq(0.5, n + 0.5, by = 1),
    col = "white",
    lwd = 1
  )
  
  axis(
    side = 1,
    at = seq_len(n),
    labels = colnames(mat),
    las = 2,
    tick = FALSE,
    line = 0.0,
    cex.axis = cex_axis
  )
  
  axis(
    side = 2,
    at = seq_len(n),
    labels = rev(rownames(mat)),
    las = 2,
    tick = FALSE,
    line = 0.0,
    cex.axis = cex_axis
  )
  
  mat_plot <- mat[n:1, , drop = FALSE]
  for (i in seq_len(nrow(mat_plot))) {
    for (j in seq_len(ncol(mat_plot))) {
      val <- mat_plot[i, j]
      txt_col <- if (abs(val) >= 0.55) "white" else "black"
      text(
        x = j,
        y = i,
        labels = sprintf("%.2f", val),
        cex = cex_cell,
        col = txt_col
      )
    }
  }
  
  box(col = "grey40", lwd = 1.2)
  
  # ----------------------------------------------------------------
  # Color key panel
  # ----------------------------------------------------------------
  par(
    mar = c(key_bottom_mar, key_left_mar, heatmap_top_mar, key_right_mar)
  )
  
  key_y <- seq(zlim[1], zlim[2], length.out = n_colors)
  key_z <- matrix(rep(key_y, each = 2), nrow = 2, byrow = TRUE)
  
  image(
    x = c(0, 1),
    y = key_y,
    z = key_z,
    col = cols,
    breaks = breaks,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  
  axis(
    side = 4,
    at = seq(-1, 1, by = 0.25),
    labels = sprintf("%.2f", seq(-1, 1, by = 0.25)),
    las = 2,
    cex.axis = key_axis_cex
  )
  
  mtext(key_label, side = 4, line = 3.8, cex = 1.00)
  box(col = "grey40", lwd = 1.0)
}

# ------------------------------------------------------------------
# Export Figure S4
# ------------------------------------------------------------------
make_figure_S4 <- function(
    tiff_file,
    pdf_file,
    cor_mat,
    label_map = pretty_labels,
    fig_width_in = 14.0,
    fig_height_in = 10.5,
    dpi = 600,
    cex_axis = 0.96,
    cex_cell = 0.82,
    debug_margins = FALSE
) {
  tiff(
    filename = tiff_file,
    width = fig_width_in,
    height = fig_height_in,
    units = "in",
    res = dpi,
    compression = "lzw"
  )
  plot_cor_heatmap_pub(
    cor_mat = cor_mat,
    label_map = label_map,
    cex_axis = cex_axis,
    cex_cell = cex_cell,
    debug_margins = debug_margins
  )
  dev.off()
  
  pdf(
    file = pdf_file,
    width = fig_width_in,
    height = fig_height_in,
    useDingbats = FALSE
  )
  plot_cor_heatmap_pub(
    cor_mat = cor_mat,
    label_map = label_map,
    cex_axis = cex_axis,
    cex_cell = cex_cell,
    debug_margins = debug_margins
  )
  dev.off()
}

make_figure_S4(
  tiff_file = file.path(OUTDIR, "Figure_S4_metrics_spearman_heatmap_publication_600dpi_autofit.tiff"),
  pdf_file  = file.path(OUTDIR, "Figure_S4_metrics_spearman_heatmap_publication_autofit.pdf"),
  cor_mat   = cor_mat,
  label_map = pretty_labels,
  fig_width_in  = 14.0,
  fig_height_in = 10.5,
  dpi = 600,
  cex_axis = 0.96,
  cex_cell = 0.82,
  debug_margins = TRUE
)
