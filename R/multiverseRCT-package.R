#' multiverseRCT – Analyse "what-if" research paths in randomised trials
#'
#' Tools to specify, run, and summarise \emph{multiverse} analyses
#' in behavioural RCTs.  The package provides:
#' \itemize{
#'   \item declarative multiverse specification with tidy syntax;
#'   \item batch execution with progress reporting and caching;
#'   \item tidy outputs that integrate with the \pkg{broom} and \pkg{dplyr} ecosystems;
#'   \item visual summaries such as effect-size rainclouds and specification curves.
#' }
#'
#' @docType package
#' @name multiverseRCT
#'
#' @import dplyr ggplot2 tidyr rlang cli tibble future.apply
#' @import gbm mgcv
#' @importFrom stats coef na.omit reorder
#' @importFrom utils globalVariables
"_PACKAGE"

## ---- globalVariables ---------------------------------------------
utils::globalVariables(c(
  ".data", ".model_name", ".prep_name", "estimate", "treatment",
  "post_test_score", "pre_test_score"
))
