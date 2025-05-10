# ------------------------------------------------------------------
#  multiverseR – Quick multiverse analyses for experimental data
# ------------------------------------------------------------------
#  Put this file in  multiverseR/R/multiverse.R
# ------------------------------------------------------------------

# ---- Imports (stand-alone block) ---------------------------------
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import cli
#' @import rlang
#' @importFrom purrr pmap_dfr
#' @importFrom future.apply future_lapply
NULL
# ------------------------------------------------------------------

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# ====================  CORE FUNCTION ==============================

#' Run a multiverse of analytic decisions
#'
#' @inheritParams defaults_prep
#' @param data        Data frame.
#' @param outcome     Bare **unquoted** outcome variable.
#' @param treatment   Bare **unquoted** treatment indicator (0/1 or factor).
#' @param covariates  Optional bare names or one-sided formula.
#' @param preprocessing Named list of preprocessing functions.
#' @param model_specs   Named list of model specifications.
#' @param id          Optional id column (bare name).
#' @param parallel    Use `future.apply` if `TRUE`.
#' @param ...         Extra args passed to custom `fit()` functions.
#'
#' @return An object of class **multiverse**.
#' @export
multiverse <- function(data,
                       outcome,
                       treatment,
                       covariates   = NULL,
                       preprocessing = defaults_prep(),
                       model_specs  = defaults_models(),
                       id           = NULL,
                       parallel     = FALSE,
                       ...) {

  stopifnot(is.list(preprocessing), is.list(model_specs))

  outcome   <- rlang::ensym(outcome)
  treatment <- rlang::ensym(treatment)
  id        <- rlang::enexpr(id)

  # -- Cartesian product of analytic choices -----------------------
  grid <- tidyr::expand_grid(
    .prep_name  = names(preprocessing),
    .model_name = names(model_specs)
  )

  run_fun <- if (parallel) future.apply::future_lapply else lapply

  grid$res <- run_fun(seq_len(nrow(grid)), function(i) {

    prep_fn   <- preprocessing[[ grid$.prep_name[i] ]]
    model_def <- model_specs [[ grid$.model_name[i] ]]

    dat_i <- prep_fn(data)

    # ----- build formula ------------------------------------------
    rhs <- switch(
      TRUE,
      is.null(covariates) ~ rlang::as_string(treatment),
      inherits(covariates, "formula") ~ paste(
        rlang::as_string(treatment),
        rlang::f_rhs(covariates), sep = " + "),
      TRUE ~ paste(
        rlang::as_string(treatment),
        paste(dplyr::as_label(rlang::ensyms(covariates)),
              collapse = " + "),
        sep = " + ")
    )
    fmla <- stats::as.formula(paste(rlang::as_string(outcome), "~", rhs))

    # ----- fit model ---------------------------------------------
    fit <- tryCatch({
      if (inherits(model_def, "formula")) {
        stats::lm(model_def, data = dat_i)
      } else if (is.list(model_def) && !is.null(model_def$fit)) {
        model_def$fit(formula = model_def$formula %||% fmla,
                      data    = dat_i, ...)
      } else stop("Bad entry in model_specs.")
    }, error = function(e) {
      warning("Model failed for ",
              grid$.prep_name[i], " × ", grid$.model_name[i], ": ",
              conditionMessage(e))
      NULL
    })

    # ----- estimate ----------------------------------------------
    est <- tryCatch({
      if (is.null(fit)) return(NA_real_)
      if (inherits(fit, "lm")) {
        coef(summary(fit))[rlang::as_string(treatment), "Estimate"]
      } else if (inherits(fit, "gbm")) {
        preds <- stats::predict(fit, dat_i, n.trees = fit$n.trees)
        mean(preds[dat_i[[rlang::as_string(treatment)]] == 1], na.rm = TRUE) -
          mean(preds[dat_i[[rlang::as_string(treatment)]] == 0], na.rm = TRUE)
      } else {
        preds <- stats::predict(fit, dat_i)
        mean(preds[dat_i[[rlang::as_string(treatment)]] == 1], na.rm = TRUE) -
          mean(preds[dat_i[[rlang::as_string(treatment)]] == 0], na.rm = TRUE)
      }
    }, error = function(e) NA_real_)

    list(fit = fit, estimate = est, n = nrow(dat_i))
  })

  class(grid) <- c("multiverse", class(grid))
  grid
}

# ====================  S3 METHODS =================================

#' @export
print.multiverse <- function(x, ...) {
  cli::cli_h1("multiverse object")
  cli::cli_ul(c(
    paste0("Universes analysed : ", nrow(x)),
    paste0("Preprocessing opts : ", length(unique(x$.prep_name))),
    paste0("Model specs        : ", length(unique(x$.model_name)))
  ))
  invisible(x)
}

#' @export
summary.multiverse <- function(object, ...) {
  purrr::pmap_dfr(
    list(object$.prep_name, object$.model_name, object$res),
    function(prep, mod, res) {
      if (!is.list(res)) res <- list()
      tibble::tibble(
        .prep_name  = prep,
        .model_name = mod,
        estimate    = res$estimate %||% NA_real_,
        n           = res$n        %||% NA_integer_
      )
    }
  )
}

#' @export
plot.multiverse <- function(x,
                            type = c("spec_curve", "distribution"),
                            ...) {
  type <- match.arg(type)
  df   <- summary(x)

  if (type == "spec_curve") {
    ggplot2::ggplot(df,
                    ggplot2::aes(x = reorder(.prep_name, estimate), y = estimate)) +
      ggplot2::geom_point(ggplot2::aes(shape = .model_name),
                          size = 2, colour = "black") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = "Preprocessing option",
                    y = "Estimated treatment effect",
                    shape = "Model",
                    title = "Specification curve")
  } else {
    ggplot2::ggplot(df,
                    ggplot2::aes(x = .model_name, y = estimate)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL, y = "Estimated treatment effect",
                    title = "Effect distribution across choices")
  }
}

# ====================  DEFAULT CHOICE SETS ========================

#' Default preprocessing pipelines
#' @export
defaults_prep <- function() {
  list(
    none             = dplyr::as_tibble,
    remove_strict    = function(df) remove_outliers(df, "post_test_score", 0.05),
    remove_lenient   = function(df) remove_outliers(df, "post_test_score", 0.01),
    winsorize        = function(df) winsorize(df, "post_test_score"),
    mean_impute      = mean_impute,
    median_impute    = median_impute,
    standardize      = standardize,
    log_transform    = transform_skewed
  )
}

#' Default model specifications (LM, RF, GAM, GBM)
#' @export
defaults_models <- function() {
  list(
    lm_base = ~ treatment,
    lm_covs = ~ treatment + pre_test_score + age,
    rf_ext  = list(
      formula = ~ post_test_score ~ treatment + pre_test_score +
        age + education_years + engagement_score + motivation_score,
      fit = function(formula, data, ...) {
        randomForest::randomForest(formula, data = data,
                                   na.action = na.omit, ...)
      }
    ),
    gam_ext = list(
      formula = ~ post_test_score ~ treatment + s(pre_test_score) +
        s(age) + s(engagement_score) +
        education_years + motivation_score,
      fit = function(formula, data, ...) mgcv::gam(formula, data = data, ...)
    ),
    gbm_ext = list(
      formula = ~ post_test_score ~ treatment + pre_test_score +
        age + education_years + engagement_score + motivation_score,
      fit = function(formula, data, n.trees = 100, ...) {
        gbm::gbm(formula, data = data, distribution = "gaussian",
                 n.trees = n.trees, verbose = FALSE, ...)
      }
    )
  )
}

# ====================  PRE-PROC HELPERS ===========================

#' Remove rows outside central proportion *p*
remove_outliers <- function(df, column, p = 0.05) {
  q <- stats::quantile(df[[column]], probs = c(p, 1 - p), na.rm = TRUE)
  dplyr::filter(df, dplyr::between(.data[[column]], q[1], q[2]))
}

#' Winsorise a numeric column
#' @export
winsorize <- function(df, column, p = 0.025) {
  q <- stats::quantile(df[[column]], probs = c(p, 1 - p), na.rm = TRUE)
  dplyr::mutate(df,
                dplyr::across(dplyr::all_of(column),
                              ~ dplyr::case_when(
                                .x > q[2] ~ q[2],
                                .x < q[1] ~ q[1],
                                TRUE      ~ .x)))
}

#' Mean imputation within treatment arms
#' @export
mean_impute <- function(df) {
  df %>%
    dplyr::group_by(treatment) %>%
    dplyr::mutate(dplyr::across(where(is.numeric),
                                ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) %>%
    dplyr::ungroup()
}

#' Median imputation within treatment arms
#' @export
median_impute <- function(df) {
  df %>%
    dplyr::group_by(treatment) %>%
    dplyr::mutate(dplyr::across(where(is.numeric),
                                ~ ifelse(is.na(.x), stats::median(.x, na.rm = TRUE), .x))) %>%
    dplyr::ungroup()
}

standardize <- function(df) {
  dplyr::mutate(df,
                dplyr::across(c(pre_test_score, post_test_score), scale))
}

transform_skewed <- function(df) {
  dplyr::mutate(df,
                dplyr::across(c(pre_test_score, post_test_score),
                              ~ scale(log1p(.x - min(.x, na.rm = TRUE) + 1))))
}

# ---- END OF FILE -------------------------------------------------
