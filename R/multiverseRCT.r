# ------------------------------------------------------------------
#  multiverseRCT  –  Multiverse analyses for behavioural RCT data
#  (ENHANCED VERSION WITH IMPROVEMENTS)
# ------------------------------------------------------------------
#  Save this file as  R/multiverseRCT.R
# ------------------------------------------------------------------

# ======= Package overview + selective imports =====================

#' multiverseRCT – explore analytic multiverses for randomised trials
#'
#' Supplies a small DSL to (a) preprocess data in multiple ways, (b) combine
#' those versions with multiple model specifications, and (c) summarise /
#' visualise the result grid.  No variable names are hard-coded: users pass
#' `outcome`, `treatment`, and optional `covariates`.
#'
#' @docType package
#' @name multiverseRCT
#'
#' @import dplyr tidyr ggplot2 cli rlang tibble
#' @importFrom purrr pmap_dfr map_dfr
#' @importFrom future.apply future_lapply
#' @importFrom stats coef predict quantile median sd binomial poisson na.omit
"_PACKAGE"

# Define global variables to avoid R CMD check NOTES
utils::globalVariables(c(".data", ".model_name", ".prep_name", "estimate",
                         "decision", "effect_size", "sort_var", "significant",
                         "lower", "upper", "p.value"))

# ------------------------------------------------------------------
#  Internal null-coalescing operator
# ------------------------------------------------------------------
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------------
#  Input validation functions
# ------------------------------------------------------------------
#' Validate inputs for multiverse analysis
#' @keywords internal
#' @noRd
validate_inputs <- function(data, treatment_sym, outcome_sym, covars_chr) {
  # Check if variables exist in data
  treat_str <- rlang::as_string(treatment_sym)
  outcome_str <- rlang::as_string(outcome_sym)

  if (!treat_str %in% names(data))
    stop("Treatment variable '", treat_str, "' not found in data")
  if (!outcome_str %in% names(data))
    stop("Outcome variable '", outcome_str, "' not found in data")

  # Check if covariates exist
  if (length(covars_chr) > 0) {
    missing_covars <- setdiff(covars_chr, names(data))
    if (length(missing_covars) > 0)
      stop("Covariates not found in data: ", paste(missing_covars, collapse = ", "))
  }

  # Check treatment is binary
  treat_var <- data[[treat_str]]
  if (is.numeric(treat_var)) {
    unique_vals <- unique(stats::na.omit(treat_var))
    if (length(unique_vals) != 2 || !all(unique_vals %in% c(0, 1)))
      warning("Treatment variable should be binary (0/1)")
  } else if (is.factor(treat_var)) {
    if (length(levels(treat_var)) != 2)
      warning("Treatment factor should have exactly 2 levels")
  } else {
    warning("Treatment variable should be binary numeric or factor")
  }

  # Check outcome type (for informational purposes)
  outcome_var <- data[[outcome_str]]
  if (is.numeric(outcome_var)) {
    if (all(outcome_var %in% c(0, 1), na.rm = TRUE)) {
      message("Outcome appears to be binary. Consider using outcome_type='binary'")
    } else if (all(outcome_var >= 0 & round(outcome_var) == outcome_var, na.rm = TRUE)) {
      message("Outcome appears to be count data. Consider using outcome_type='count'")
    }
  }

  invisible(data)
}

# ------------------------------------------------------------------
#  Effect extraction functions
# ------------------------------------------------------------------
#' Extract treatment effect from various model types
#' @keywords internal
#' @noRd
extract_effect <- function(fit, treatment, data) {
  if (is.null(fit)) return(NA_real_)

  tryCatch({
    if (inherits(fit, "lm") || inherits(fit, "glm")) {
      # For linear/generalized linear models, extract coefficient
      return(stats::coef(fit)[treatment])
    } else if (inherits(fit, "randomForest") || inherits(fit, "ranger")) {
      # For random forest models, use prediction differences
      # First, create copies of the data with treatment = 0/1
      data_0 <- data
      data_1 <- data
      data_0[[treatment]] <- 0
      data_1[[treatment]] <- 1

      # Predict outcomes under both scenarios
      pred_0 <- predict(fit, data_0)
      pred_1 <- predict(fit, data_1)

      # Return average treatment effect
      return(mean(pred_1 - pred_0, na.rm = TRUE))
    } else {
      # For other model types, use generic mean difference in predictions
      preds <- predict(fit, data)
      treat_var <- data[[treatment]]
      m1 <- mean(preds[treat_var == 1], na.rm = TRUE)
      m0 <- mean(preds[treat_var == 0], na.rm = TRUE)
      return(m1 - m0)
    }
  }, error = function(e) {
    warning("Failed to extract effect: ", e$message)
    return(NA_real_)
  })
}

# ====================  CORE FUNCTION ==============================

#' Run a multiverse of analytic decisions
#'
#' @param data Data frame.
#' @param outcome Bare outcome variable.
#' @param treatment Bare treatment indicator (0/1 or factor).
#' @param covariates Optional character vector or one-sided formula
#'   with covariates, e.g. `~ age + baseline`.
#' @param preprocessing Named list of preprocessing functions.
#' @param model_specs Named list of model definitions.  If `NULL`, models
#'   are created automatically based on `outcome_type`.
#' @param outcome_type Type of outcome variable: "continuous", "binary", or "count".
#' @param parallel Logical, use `future_lapply()` if `TRUE`. Can also be a number
#'   specifying the number of worker cores to use.
#' @param extract_fun Custom function to extract treatment effects from models.
#' @param ... Extra args passed to custom `fit()` functions.
#'
#' @return A tibble of class `"multiverse_rct"`.
#' @export
#'
#' @examples
#' # Create a simple example dataset
#' set.seed(123)
#' n <- 30
#' example_data <- data.frame(
#'   treatment = rep(c(0, 1), each = n/2),
#'   covariate = rnorm(n),
#'   outcome = rnorm(n)
#' )
#' example_data$outcome <- example_data$outcome + 0.5 * example_data$treatment
#'
#' # Basic usage with continuous outcome
#' mv <- multiverse_rct(
#'   data = example_data,
#'   outcome = outcome,
#'   treatment = treatment,
#'   covariates = ~ covariate
#' )
#'
#' # Summarize results
#' summary(mv)
#'
#' \dontrun{
#' # Examples using the included rct_data dataset
#'
#' # With continuous outcome
#' mv <- multiverse_rct(
#'   data = rct_data,
#'   outcome = outcome_score,
#'   treatment = treatment_group,
#'   covariates = ~ baseline + age
#' )
#'
#' # With binary outcome
#' mv_binary <- multiverse_rct(
#'   data = rct_data,
#'   outcome = success,
#'   treatment = treatment_group,
#'   outcome_type = "binary"
#' )
#'
#' # With custom preprocessing
#' my_preprocessing <- defaults_prep()
#' my_preprocessing$log_transform <- function(df) {
#'   df$outcome_score <- log(df$outcome_score + 1)
#'   return(df)
#' }
#'
#' mv_custom <- multiverse_rct(
#'   data = rct_data,
#'   outcome = outcome_score,
#'   treatment = treatment_group,
#'   preprocessing = my_preprocessing
#' )
#' }
multiverse_rct <- function(data,
                           outcome,
                           treatment,
                           covariates = NULL,
                           preprocessing = defaults_prep(),
                           model_specs = NULL,
                           outcome_type = c("continuous", "binary", "count"),
                           parallel = FALSE,
                           extract_fun = NULL,
                           ...) {

  outcome_sym <- rlang::ensym(outcome)
  treatment_sym <- rlang::ensym(treatment)
  outcome_type <- match.arg(outcome_type)

  # Convert covariates to character vector
  covars_chr <- character()
  if (inherits(covariates, "formula")) {
    covars_chr <- all.vars(rlang::f_rhs(covariates))
  } else if (!is.null(covariates)) {
    covars_chr <- as.character(covariates)
  }

  # Validate inputs
  validate_inputs(data, treatment_sym, outcome_sym, covars_chr)

  # Default models if none supplied, based on outcome type
  if (is.null(model_specs)) {
    model_specs <- defaults_models_by_type(
      outcome = rlang::as_string(outcome_sym),
      treatment = rlang::as_string(treatment_sym),
      covars = covars_chr,
      outcome_type = outcome_type
    )
  }

  # Use custom effect extraction function if provided
  if (is.null(extract_fun)) {
    extract_fun <- function(fit, dat) {
      extract_effect(
        fit = fit,
        treatment = rlang::as_string(treatment_sym),
        data = dat
      )
    }
  }

  # Set up parallel processing if requested
  if (isTRUE(parallel) || (is.numeric(parallel) && parallel > 1)) {
    if (!requireNamespace("future", quietly = TRUE)) {
      warning("The 'future' package is required for parallel processing. Using sequential processing instead.")
      run_fun <- lapply
    } else {
      # Set up a default plan if none exists or if cores specified
      if (is.numeric(parallel)) {
        future::plan(future::multisession, workers = parallel)
      } else if (is.null(future::plan())) {
        future::plan(future::multisession,
                     workers = min(parallel::detectCores() - 1, 4))
      }
      run_fun <- future.apply::future_lapply
    }
  } else {
    run_fun <- lapply
  }

  # Create analysis grid
  grid <- tidyr::expand_grid(
    .prep_name = names(preprocessing),
    .model_name = names(model_specs)
  )

  # Run all models in the multiverse
  grid$res <- run_fun(seq_len(nrow(grid)), function(i) {

    dat_i <- preprocessing[[ grid$.prep_name[i] ]](data)
    spec_i <- model_specs[[ grid$.model_name[i] ]]

    # Build formula
    if (inherits(spec_i, "formula")) {
      fmla <- spec_i
      fitfun <- stats::lm
    } else {
      fmla <- spec_i$build(
        outcome = rlang::as_string(outcome_sym),
        treatment = rlang::as_string(treatment_sym),
        covars = covars_chr
      )
      fitfun <- spec_i$fit %||% stats::lm
    }

    # Fit model
    fit <- tryCatch(
      fitfun(fmla, data = dat_i, ...),
      error = function(e) {
        warning(paste("Error in", grid$.prep_name[i], "-",
                      grid$.model_name[i], ":", e$message))
        NULL
      }
    )

    # Extract effect estimate
    est <- if (is.null(fit)) NA_real_ else extract_fun(fit, dat_i)

    # Get standard error if available
    se <- tryCatch({
      if (inherits(fit, "lm") || inherits(fit, "glm")) {
        summary(fit)$coefficients[rlang::as_string(treatment_sym), "Std. Error"]
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)

    # Get p-value if available
    p_value <- tryCatch({
      if (inherits(fit, "lm") || inherits(fit, "glm")) {
        summary(fit)$coefficients[rlang::as_string(treatment_sym), "Pr(>|t|)"]
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)

    list(
      fit = fit,
      estimate = est,
      std.error = se,
      p.value = p_value,
      n = nrow(dat_i),
      dat = dat_i  # Store the preprocessed data for later use
    )
  })

  result <- grid
  attr(result, "outcome") <- rlang::as_string(outcome_sym)
  attr(result, "treatment") <- rlang::as_string(treatment_sym)
  attr(result, "covariates") <- covars_chr
  attr(result, "outcome_type") <- outcome_type
  class(result) <- c("multiverse_rct", class(grid))
  result
}

# ====================  S3 METHODS =================================

#' @export
print.multiverse_rct <- function(x, ...) {
  cli::cli_h1("multiverseRCT object")

  # Extract summary information
  n_universes <- nrow(x)
  n_preps <- length(unique(x$.prep_name))
  n_models <- length(unique(x$.model_name))

  # Get outcome and treatment variable names
  outcome <- attr(x, "outcome") %||% "unknown"
  treatment <- attr(x, "treatment") %||% "unknown"
  outcome_type <- attr(x, "outcome_type") %||% "continuous"

  # Calculate summary stats
  summ <- summary(x)
  mean_effect <- mean(summ$estimate, na.rm = TRUE)
  median_effect <- stats::median(summ$estimate, na.rm = TRUE)
  sd_effect <- stats::sd(summ$estimate, na.rm = TRUE)
  n_positive <- sum(summ$estimate > 0, na.rm = TRUE)
  n_negative <- sum(summ$estimate < 0, na.rm = TRUE)
  n_significant <- sum(summ$p.value < 0.05, na.rm = TRUE)

  # Print information
  cli::cli_ul(c(
    paste0("Outcome variable      : ", outcome, " (", outcome_type, ")"),
    paste0("Treatment variable    : ", treatment),
    paste0("Universes analyzed    : ", n_universes),
    paste0("Pre-processing options: ", n_preps),
    paste0("Model specifications  : ", n_models),
    paste0("Mean effect estimate  : ", round(mean_effect, 4)),
    paste0("Median effect estimate: ", round(median_effect, 4)),
    paste0("Effect SD (vibration) : ", round(sd_effect, 4)),
    paste0("Direction consistency : ",
           round(max(n_positive, n_negative) / (n_positive + n_negative) * 100, 1), "%"),
    paste0("Significant results   : ", n_significant, "/", n_universes,
           " (", round(n_significant/n_universes*100, 1), "%)")
  ))

  invisible(x)
}

#' @export
summary.multiverse_rct <- function(object, ...) {
  purrr::map_dfr(seq_len(nrow(object)), function(i) {
    prep <- object$.prep_name[i]
    mod <- object$.model_name[i]
    res <- object$res[[i]]

    if (!is.list(res)) res <- list()

    tibble::tibble(
      .prep_name = prep,
      .model_name = mod,
      estimate = res$estimate %||% NA_real_,
      std.error = res$std.error %||% NA_real_,
      p.value = res$p.value %||% NA_real_,
      n = res$n %||% NA_integer_,
      converged = !is.null(res$fit)
    )
  })
}

#' @export
plot.multiverse_rct <- function(x,
                                type = c("spec_curve", "distribution",
                                         "heatmap", "forest", "significance"),
                                signif_level = 0.05,
                                sort_by = c("estimate", "p.value", "preprocessing", "model"),
                                ...) {
  type <- match.arg(type)
  sort_by <- match.arg(sort_by)
  df <- summary(x)

  # Add formatted decision labels
  df$decision <- paste(df$.prep_name, "-", df$.model_name)

  # Add significance indicator
  df$significant <- df$p.value < signif_level

  # Determine sort order
  if (sort_by == "estimate") {
    df$sort_var <- df$estimate
  } else if (sort_by == "p.value") {
    df$sort_var <- df$p.value
  } else if (sort_by == "preprocessing") {
    df$sort_var <- paste0(df$.prep_name, "_", df$.model_name)
  } else {
    df$sort_var <- paste0(df$.model_name, "_", df$.prep_name)
  }

  if (type == "spec_curve") {
    # Specification curve with optional confidence intervals
    p <- ggplot2::ggplot(df,
                         ggplot2::aes(x = reorder(decision, sort_var),
                                      y = estimate)) +
      ggplot2::geom_point(ggplot2::aes(shape = .model_name,
                                       color = significant)) +
      ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = "Analytic Decision",
                    y = "Estimated Treatment Effect",
                    shape = "Model",
                    color = "p < 0.05",
                    title = "Specification Curve") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed")

    # Add confidence intervals if standard errors available
    if (!all(is.na(df$std.error))) {
      df$lower <- df$estimate - 1.96 * df$std.error
      df$upper <- df$estimate + 1.96 * df$std.error
      p <- p + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                                       height = 0.2, alpha = 0.5)
    }

    return(p)

  } else if (type == "distribution") {
    # Distribution of effects by model type
    ggplot2::ggplot(df,
                    ggplot2::aes(x = .model_name, y = estimate, fill = .model_name)) +
      ggplot2::geom_violin(alpha = 0.7) +
      ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
      ggplot2::geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = "Model Specification",
                    y = "Estimated Treatment Effect",
                    title = "Effect Distribution Across Model Choices") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed")

  } else if (type == "heatmap") {
    # Heatmap of effects
    ggplot2::ggplot(df,
                    ggplot2::aes(x = .model_name, y = .prep_name, fill = estimate)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                    midpoint = 0) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "Treatment Effect Heatmap",
                    x = "Model Specification",
                    y = "Pre-processing Method",
                    fill = "Effect Size")

  } else if (type == "forest") {
    # Forest plot of all effects
    ggplot2::ggplot(df,
                    ggplot2::aes(x = estimate,
                                 y = reorder(decision, sort_var),
                                 color = significant)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::theme_minimal() +
      ggplot2::labs(y = NULL,
                    x = "Treatment Effect",
                    color = "p < 0.05",
                    title = "Forest Plot of Treatment Effects")

  } else if (type == "significance") {
    # Plot focusing on significance patterns
    ggplot2::ggplot(df,
                    ggplot2::aes(x = estimate, y = -log10(p.value),
                                 color = significant)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      ggplot2::geom_hline(yintercept = -log10(signif_level),
                          linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::annotate("text", x = min(df$estimate, na.rm = TRUE),
                        y = -log10(signif_level) + 0.2,
                        label = paste("p =", signif_level),
                        hjust = 0) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Effect Estimate",
                    y = "-log10(p-value)",
                    color = paste("p <", signif_level),
                    title = "Significance Plot")
  }
}

# ====================  ADDITIONAL FUNCTIONS =======================

#' Filter multiverse results
#'
#' Select a subset of the multiverse based on various criteria.
#'
#' @param x A multiverse_rct object
#' @param min_effect Minimum effect size to include
#' @param max_effect Maximum effect size to include
#' @param min_n Minimum sample size to include
#' @param prep Preprocessing method(s) to include
#' @param model Model specification(s) to include
#' @param significant_only If TRUE, only include models with p < 0.05
#'
#' @return A filtered multiverse_rct object
#' @export
filter_multiverse <- function(x, min_effect = NULL, max_effect = NULL,
                              min_n = NULL, prep = NULL, model = NULL,
                              significant_only = FALSE) {
  df <- summary(x)

  if (!is.null(min_effect))
    df <- df[df$estimate >= min_effect, ]
  if (!is.null(max_effect))
    df <- df[df$estimate <= max_effect, ]
  if (!is.null(min_n))
    df <- df[df$n >= min_n, ]
  if (!is.null(prep))
    df <- df[df$.prep_name %in% prep, ]
  if (!is.null(model))
    df <- df[df$.model_name %in% model, ]
  if (significant_only)
    df <- df[df$p.value < 0.05, ]

  # Subset original multiverse object
  idx <- paste(df$.prep_name, df$.model_name)
  full_idx <- paste(x$.prep_name, x$.model_name)

  result <- x[full_idx %in% idx, ]
  attributes(result) <- attributes(x)
  class(result) <- class(x)
  return(result)
}

#' Extract models from multiverse results
#'
#' Get the fitted model objects from specific specifications.
#'
#' @param x A multiverse_rct object
#' @param prep_name Preprocessing method name(s) or NULL for all
#' @param model_name Model specification name(s) or NULL for all
#'
#' @return A list of fitted model objects
#' @export
extract_models <- function(x, prep_name = NULL, model_name = NULL) {
  indices <- rep(TRUE, nrow(x))

  if (!is.null(prep_name))
    indices <- indices & (x$.prep_name %in% prep_name)
  if (!is.null(model_name))
    indices <- indices & (x$.model_name %in% model_name)

  models <- lapply(x$res[indices], function(res) res$fit)
  names(models) <- paste(x$.prep_name[indices], x$.model_name[indices], sep = "-")
  models
}

#' Sensitivity analysis for multiverse results
#'
#' Assess the stability of findings across the analytic multiverse.
#'
#' @param x A multiverse_rct object
#'
#' @return A list of sensitivity metrics
#' @export
sensitivity_analysis <- function(x) {
  df <- summary(x)

  # Calculate overall summary
  overall_mean <- mean(df$estimate, na.rm = TRUE)
  overall_sd <- stats::sd(df$estimate, na.rm = TRUE)
  overall_median <- stats::median(df$estimate, na.rm = TRUE)

  # Effect by preprocessing method
  prep_effects <- tapply(df$estimate, df$.prep_name,
                         function(x) c(mean = mean(x, na.rm = TRUE),
                                       sd = stats::sd(x, na.rm = TRUE)))

  # Effect by model type
  model_effects <- tapply(df$estimate, df$.model_name,
                          function(x) c(mean = mean(x, na.rm = TRUE),
                                        sd = stats::sd(x, na.rm = TRUE)))

  # Significance consistency
  sig_consistency <- mean(df$p.value < 0.05, na.rm = TRUE)

  # Sign consistency (proportion of results with same sign)
  sign_consistency <- mean(sign(df$estimate) == sign(overall_median), na.rm = TRUE)

  # Vibration of effects ratio
  voe <- overall_sd / abs(overall_median)

  # Prepare results
  result <- list(
    overall = list(
      mean = overall_mean,
      median = overall_median,
      sd = overall_sd
    ),
    by_preprocessing = prep_effects,
    by_model = model_effects,
    sign_consistency = sign_consistency,
    significance_rate = sig_consistency,
    vibration_of_effects = voe,
    interpretation = if (voe < 0.5) {
      "Findings appear robust to analytic decisions"
    } else if (voe < 1) {
      "Moderate sensitivity to analytic decisions"
    } else {
      "High sensitivity to analytic decisions"
    }
  )

  class(result) <- "multiverse_sensitivity"
  return(result)
}

#' @export
print.multiverse_sensitivity <- function(x, ...) {
  cli::cli_h1("Multiverse Sensitivity Analysis")

  # Print overall effects
  cli::cli_h2("Overall Effect Estimates")
  cli::cli_ul(c(
    paste0("Mean estimate   : ", round(x$overall$mean, 4)),
    paste0("Median estimate : ", round(x$overall$median, 4)),
    paste0("Standard dev.   : ", round(x$overall$sd, 4))
  ))

  # Print consistency metrics
  cli::cli_h2("Consistency Metrics")
  cli::cli_ul(c(
    paste0("Sign consistency     : ", round(x$sign_consistency * 100, 1), "%"),
    paste0("Significance rate    : ", round(x$significance_rate * 100, 1), "%"),
    paste0("Vibration of effects : ", round(x$vibration_of_effects, 3)),
    paste0("Interpretation       : ", x$interpretation)
  ))

  # Print effects by preprocessing
  cli::cli_h2("Effects by Preprocessing Method")
  for (i in seq_along(x$by_preprocessing)) {
    method <- names(x$by_preprocessing)[i]
    stats <- x$by_preprocessing[[i]]
    cli::cli_text(paste0(method, ": mean = ", round(stats["mean"], 4),
                         ", SD = ", round(stats["sd"], 4)))
  }

  # Print effects by model
  cli::cli_h2("Effects by Model Specification")
  for (i in seq_along(x$by_model)) {
    model <- names(x$by_model)[i]
    stats <- x$by_model[[i]]
    cli::cli_text(paste0(model, ": mean = ", round(stats["mean"], 4),
                         ", SD = ", round(stats["sd"], 4)))
  }

  invisible(x)
}

# ====================  DEFAULTS (ENHANCED) ========================

#' Default preprocessing pipelines
#'
#' @param na_threshold Maximum proportion of missing values allowed before
#'   excluding a variable (used in the `drop_sparse` method)
#'
#' @return A named list of preprocessing functions
#' @export
defaults_prep <- function(na_threshold = 0.2) {
  list(
    none = identity,

    winsorize = winsorize,

    mean_impute = mean_impute,

    median_impute = median_impute,

    # Added: complete case analysis
    complete_cases = function(df) {
      stats::na.omit(df)
    },

    # Added: drop sparse columns then complete cases
    drop_sparse_then_complete = function(df) {
      # Identify columns with too many NAs
      na_prop <- sapply(df, function(x) mean(is.na(x)))
      keep_cols <- na_prop <= na_threshold

      # Only try to keep demographic columns that actually exist in the data
      demo_cols <- c("treatment", "arm", "group", "id", "gender", "sex", "age")
      existing_demo_cols <- intersect(demo_cols, names(df))
      if (length(existing_demo_cols) > 0) {
        keep_cols[existing_demo_cols] <- TRUE
      }

      # Subset and do complete case analysis
      df_subset <- df[, keep_cols, drop = FALSE]
      stats::na.omit(df_subset)
    },

    # Added: standardize numeric variables
    standardize = function(df) {
      numeric_cols <- sapply(df, is.numeric)
      df[numeric_cols] <- lapply(df[numeric_cols], function(x) {
        if (length(unique(stats::na.omit(x))) > 1) {
          (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
        } else {
          x
        }
      })
      df
    },

    # Combination methods
    winsorize_and_impute = function(df) {
      df <- winsorize(df)
      mean_impute(df)
    }
  )
}

#' Default model specifications by outcome type
#'
#' @param outcome Outcome variable name
#' @param treatment Treatment variable name
#' @param covars Character vector of covariate names
#' @param outcome_type Type of outcome: "continuous", "binary", or "count"
#'
#' @return A named list of model specifications
#' @export
defaults_models_by_type <- function(outcome, treatment, covars = character(),
                                    outcome_type = c("continuous", "binary", "count")) {
  outcome_type <- match.arg(outcome_type)

  # Formula building helper
  build_formula <- function(outcome, treatment, covars = character(),
                            interaction = NULL) {
    terms <- c(treatment)

    # Add covariates
    if (length(covars) > 0) {
      # If interaction specified, add that first
      if (!is.null(interaction) && interaction %in% covars) {
        terms <- c(terms, paste0(treatment, ":", interaction))
        terms <- c(terms, setdiff(covars, interaction))
      } else {
        terms <- c(terms, covars)
      }
    }

    # Create formula string
    rhs <- paste(terms, collapse = " + ")
    stats::as.formula(paste(outcome, "~", rhs))
  }

  # Return appropriate model specifications based on outcome type
  if (outcome_type == "continuous") {
    list(
      # Basic linear model
      lm_simple = list(
        build = function(outcome, treatment, covars) {
          build_formula(outcome, treatment)
        },
        fit = stats::lm
      ),

      # Covariate-adjusted linear model
      lm_adjusted = list(
        build = build_formula,
        fit = stats::lm
      ),

      # Robust regression (if MASS is available)
      robust_reg = list(
        build = build_formula,
        fit = function(formula, data, ...) {
          if (requireNamespace("MASS", quietly = TRUE)) {
            MASS::rlm(formula, data = data, ...)
          } else {
            warning("MASS package not available, using lm instead")
            stats::lm(formula, data = data, ...)
          }
        }
      )
    )
  } else if (outcome_type == "binary") {
    list(
      # Logistic regression without covariates
      logistic_simple = list(
        build = function(outcome, treatment, covars) {
          build_formula(outcome, treatment)
        },
        fit = function(formula, data, ...) {
          stats::glm(formula, family = stats::binomial(), data = data, ...)
        }
      ),

      # Logistic regression with covariates
      logistic_adjusted = list(
        build = build_formula,
        fit = function(formula, data, ...) {
          stats::glm(formula, family = stats::binomial(), data = data, ...)
        }
      )
    )
  } else if (outcome_type == "count") {
    list(
      # Poisson regression without covariates
      poisson_simple = list(
        build = function(outcome, treatment, covars) {
          build_formula(outcome, treatment)
        },
        fit = function(formula, data, ...) {
          stats::glm(formula, family = stats::poisson(), data = data, ...)
        }
      ),

      # Poisson regression with covariates
      poisson_adjusted = list(
        build = build_formula,
        fit = function(formula, data, ...) {
          stats::glm(formula, family = stats::poisson(), data = data, ...)
        }
      ),

      # Negative binomial regression (if MASS is available)
      negbin = list(
        build = build_formula,
        fit = function(formula, data, ...) {
          if (requireNamespace("MASS", quietly = TRUE)) {
            MASS::glm.nb(formula, data = data, ...)
          } else {
            warning("MASS package not available, using poisson instead")
            stats::glm(formula, family = stats::poisson(), data = data, ...)
          }
        }
      )
    )
  }
}

#' @export
#' @rdname defaults_models_by_type
defaults_models <- function(outcome, treatment, covars = character()) {
  defaults_models_by_type(outcome, treatment, covars, "continuous")
}

# ====================  PRE-PROCESSING HELPERS (ENHANCED) ==========

#' Winsorise numeric columns
#'
#' Values below/above the (`p`, `1 – p`) quantiles are replaced by those
#' quantiles.  If `cols` is `NULL` every numeric column is winsorised.
#'
#' @param df   Data frame.
#' @param cols Character vector of columns or `NULL`.
#' @param p    Two-sided tail proportion (default 0.025).
#' @export
winsorize <- function(df, cols = NULL, p = 0.025) {
  if (is.null(cols))
    cols <- names(df)[sapply(df, is.numeric)]

  dplyr::mutate(df,
                dplyr::across(dplyr::all_of(cols), \(x) {
                  if (length(unique(stats::na.omit(x))) <= 1) return(x)
                  q <- stats::quantile(x, c(p, 1 - p), na.rm = TRUE)
                  pmin(pmax(x, q[1]), q[2])
                }))
}

#' Mean imputation for numeric columns
#' @param df Data frame
#' @param cols Character vector of columns or NULL for all numeric columns
#' @export
mean_impute <- function(df, cols = NULL) {
  if (is.null(cols))
    cols <- names(df)[sapply(df, is.numeric)]

  dplyr::mutate(df,
                dplyr::across(dplyr::all_of(cols),
                              \(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))
}

#' Median imputation for numeric columns
#' @param df Data frame
#' @param cols Character vector of columns or NULL for all numeric columns
#' @export
median_impute <- function(df, cols = NULL) {
  if (is.null(cols))
    cols <- names(df)[sapply(df, is.numeric)]

  dplyr::mutate(df,
                dplyr::across(dplyr::all_of(cols),
                              \(x) ifelse(is.na(x), stats::median(x, na.rm = TRUE), x)))
}

#' Mode imputation for categorical columns
#' @param df Data frame
#' @param cols Character vector of columns or NULL for all factor/character columns
#' @export
mode_impute <- function(df, cols = NULL) {
  if (is.null(cols))
    cols <- names(df)[sapply(df, function(x) is.factor(x) || is.character(x))]

  # Helper to find mode
  find_mode <- function(x) {
    tab <- table(x)
    names(tab)[which.max(tab)]
  }

  for (col in cols) {
    if (any(is.na(df[[col]]))) {
      mode_val <- find_mode(df[[col]])
      df[[col]][is.na(df[[col]])] <- mode_val
    }
  }

  df
}

#' Outlier detection and removal
#' @param df Data frame
#' @param cols Character vector of columns or NULL for all numeric columns
#' @param method Method for outlier detection: "z_score", "iqr", or "percentile"
#' @param threshold Threshold for outlier detection
#' @export
remove_outliers <- function(df, cols = NULL,
                            method = c("z_score", "iqr", "percentile"),
                            threshold = 3) {
  method <- match.arg(method)

  if (is.null(cols))
    cols <- names(df)[sapply(df, is.numeric)]

  # Make a copy to avoid modifying the original
  df_result <- df
  rows_to_keep <- rep(TRUE, nrow(df))

  for (col in cols) {
    x <- df[[col]]

    # Skip if not enough values
    if (sum(!is.na(x)) <= 5) next

    outlier_indices <- if (method == "z_score") {
      # Z-score method
      z <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
      abs(z) > threshold
    } else if (method == "iqr") {
      # IQR method
      q <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
      iqr <- q[2] - q[1]
      x < (q[1] - threshold * iqr) | x > (q[2] + threshold * iqr)
    } else {
      # Percentile method
      q <- stats::quantile(x, c(threshold/100, 1 - threshold/100), na.rm = TRUE)
      x < q[1] | x > q[2]
    }

    # Combine with NAs
    outlier_indices[is.na(outlier_indices)] <- FALSE

    # Update row filter
    rows_to_keep <- rows_to_keep & !outlier_indices
  }

  # Return filtered dataset
  df_result[rows_to_keep, , drop = FALSE]
}

# ---- END OF FILE -------------------------------------------------
