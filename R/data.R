#' Sample RCT dataset for examples
#'
#' A simulated randomized controlled trial dataset with treatment,
#' outcome, and covariates for demonstrating multiverseRCT functions.
#'
#' @format A data frame with 200 rows and 6 variables:
#' \describe{
#'   \item{id}{Participant ID}
#'   \item{treatment_group}{Treatment indicator (0=control, 1=treatment)}
#'   \item{baseline}{Baseline measurement}
#'   \item{age}{Participant age}
#'   \item{outcome_score}{Primary outcome variable (continuous)}
#'   \item{success}{Binary outcome (0/1)}
#'   \item{count_outcome}{Count outcome variable}
#' }
#' @source Generated for multiverseRCT examples
"rct_data"
