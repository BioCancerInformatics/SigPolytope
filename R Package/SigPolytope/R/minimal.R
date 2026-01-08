# Função mínima para teste
#' Hello World function
#' @export
hello_geocircuitry <- function() {
  cat("Hello from GeoCircuitry!\n")
}

#' Test PCA function
#' @export
test_pca <- function(x) {
  stats::prcomp(x, scale = TRUE)
}
