---
bibliography: references.bib
---
# multiverseRCT <img src="man/figures/logo.png" align="left" height="160" style="margin-right: 20px;" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/gav888/multiverseRCT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gav888/multiverseRCT/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

**multiverseRCT** is an R package for conducting multiverse analyses on randomized controlled trial (RCT) data. Instead of making a single set of analytical choices, multiverse analysis explores how different preprocessing decisions and model specifications affect research conclusions.

The package lets you:

* Apply multiple preprocessing methods (e.g., different imputation strategies, outlier handling)
* Fit multiple model specifications (adjusted/unadjusted, different model families)
* Visualize the variation in treatment effects across the multiverse
* Quantify the sensitivity of conclusions to analytic decisions

## Installation

You can install the development version of multiverseRCT like so:

``` r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("yourusername/multiverseRCT")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(multiverseRCT)

# Run a multiverse analysis with the default preprocessing and model options
mv <- multiverse_rct(
  data = your_data,
  outcome = primary_outcome,
  treatment = treatment_group,
  covariates = ~ baseline + age + stratum
)

# View summary of all results
summary(mv)

# Create a specification curve plot
plot(mv, type = "spec_curve")

# Assess sensitivity of findings
sensitivity_analysis(mv)
```

Advanced Usage

``` r
# Create a custom preprocessing pipeline
my_preprocessing <- defaults_prep()
my_preprocessing$log_transform <- function(df) {
  df$outcome_var <- log(df$outcome_var + 1)
  return(df)
}

# Run analysis with custom preprocessing
mv_custom <- multiverse_rct(
  data = your_data,
  outcome = outcome_var,
  treatment = treatment_group,
  preprocessing = my_preprocessing,
  outcome_type = "continuous"
)

# Filter to significant results only
significant_subset <- filter_multiverse(mv_custom, significant_only = TRUE)

# Compare effect distributions across models
plot(mv_custom, type = "distribution")

# Create heatmap of effects
plot(mv_custom, type = "heatmap")
```

## Key Features

-   **Variable-agnostic design**: Works with any outcome, treatment, and covariates

-   **Flexible preprocessing**: Multiple imputation strategies, outlier handling, standardization

-   **Model flexibility**: Supports continuous, binary, and count outcomes with appropriate models

-   **Visualization tools**: Specification curves, forest plots, heatmaps, and significance plots

-   **Sensitivity metrics**: Quantify the robustness of findings via vibration of effects

-   **Parallel processing**: Option for faster computation of large multiverses

## Background and Resources

Multiverse analysis addresses the issue of analytical flexibility in research by exploring how findings depend on reasonable analysis choices. Rather than choosing a single path through the garden of forking paths, multiverse analysis examines many paths simultaneously.\
\

## Licence

This package is licensed under the MIT License. See the LICENSE file for details.

## Citation

If you use multiverseRCT in your research, please cite it as:

```         
Veltri, G. A. (2025). multiverseRCT: A Toolkit for Multiverse Analyses in Randomised Controlled Trials. R package version 0.0.1.
```

## Contributing

Contributions to multiverseRCT are welcome! Please feel free to submit issues or pull requests.
