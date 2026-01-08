#' @keywords internal
.format_latent_dims_text <- function(dims, title = "Latent dimensions") {
  if (is.null(dims) || length(dims) == 0) return(NULL)
  
  lines <- paste0(seq_along(dims), " = ", dims)
  paste0(
    title, ":\n",
    paste(lines, collapse = "\n")
  )
}
