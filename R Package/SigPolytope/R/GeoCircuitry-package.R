# R/GeoCircuitry-package.R

#' @keywords internal
"_PACKAGE"

# ==============================================================================
# GeoCircuitry: Geometric Analysis of Circuitry Signatures
# ==============================================================================

#' @description
#' Provides tools for building, analyzing and visualizing geometric
#' representations of biological circuitry signatures and their regulatory
#' interactions in 3D PCA space.
#' 
#' @details
#' The main functions of the package include:
#' \itemize{
#'   \item \code{\link{build_signature_geometry}} - Build signature geometry
#'   \item \code{\link{build_regulator_geometry}} - Build regulator geometry  
#'   \item \code{\link{compute_circuitry_convergence}} - Compute convergence metrics
#'   \item \code{\link{plot_circuitry_3d}} - Interactive 3D visualization
#' }
#'
#' @author Seu Nome <seu@email.com>
#' @references 
#' Add relevant references here
#' 
#' @docType package
#' @name GeoCircuitry
#' 
#' @importFrom dplyr case_when bind_cols
#' @importFrom tibble tibble as_tibble
#' @importFrom geometry convhulln
#' @importFrom stats prcomp sd median quantile
#' @importFrom parallel makeCluster clusterEvalQ clusterExport parLapply stopCluster detectCores
#' @importFrom plotly plot_ly add_trace add_markers layout
#' 
NULL

# ==============================================================================
# Package startup message
# ==============================================================================

.onAttach <- function(libname, pkgname) {
  version <- tryCatch(
    utils::packageVersion("GeoCircuitry"),
    error = function(e) "0.1.0"
  )
  
  packageStartupMessage(
    "========================================\n",
    "GeoCircuitry ", version, "\n",
    "Geometric analysis of circuitry signatures\n",
    "========================================\n"
  )
}

# ==============================================================================
# Package unload (cleanup if needed)
# ==============================================================================

.onUnload <- function(libpath) {
  # Clean up any resources if needed
  invisible()
}

# ==============================================================================
# Global variables (if any)
# ==============================================================================

utils::globalVariables(c(
  # Variables used in dplyr pipes that might cause R CMD check warnings
  "PC1", "PC2", "PC3", 
  "barycenter_distance", "distance_implication",
  "sig_hull_vol", "int_hull_vol", "vol_ratio", "vol_implication"
))