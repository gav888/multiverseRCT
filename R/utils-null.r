# R/utils-null.R
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

