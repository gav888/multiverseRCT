# R/utils.R
# ------------------------------------------------------------------
#  Internal infix operator  %||%
#  (returns y when x is NULL, otherwise x)
# ------------------------------------------------------------------

#' @noRd
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
