################################################################################
# File: estimate.R                                                             #
# Project: crcsim                                                              #
# Created Date: 2026-08-26                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-27                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' @title Estimate the N_pop size
#'
#' @description Estimate capture probability with or without covariates.
#'
#' @details This is the controlling estimation function that determines the
#' estimation approach and coordinates data processing based on the data's
#' format. First, it validates that all the capture columns are binary (0 or 1).
#' Secondly, it validates the number of cross folds to ensure there is
#' sufficient data. We consider "sufficient" to be 50 samples per fold.
#'
#' @param data The dataset containing capture columns and optionally covariates.
#' @param n_lists The number of capture lists.
#' @param method The estimation method to be used (default is "doubly_robust").
#' @param func The function used for nuisance parameter estimation (default is
#' "logit").
#' @param nfolds The number of cross-validation folds (default is 5).
#' @param margin The margin parameter for estimation (default is 0.005).
#' @param seed The random seed for reproducibility (default is NULL).
#' @param diagnostics Controls GLM fit diagnostics: "quiet" (default),
#' "warning", or "verbose".
#' @param ... Additional arguments passed to the nuisance estimation function.
#' @return A list containing the estimation results and the total population
#' size.
#'
#' @export
estimate_n <- function(
    data,
    n_capture_cols,
    estimator_variant = log_linear_estimator(),
    seed = NULL
) {
    # validate format first: check column names/order and no NA values
    stopifnot(validate_binary_cols(data, n_capture_cols))
    if (nrow(data) < 1) {
        stop("Data must have at least one row.")
    }

    if (inherits(estimator_variant, "doubly_robust_estimator")) {
        if ((ncol(data) - n_capture_cols) < 1) {
            warning("Doubly robust estimator chosen, but no covariates found.")
        }
        return(estimate_doubly_robust(data, n_capture_cols, estimator_variant))
    } else if (inherits(estimator_variant, "log_linear_estimator")) {
        return(estimate_log_linear(data, n_capture_cols, estimator_variant))
    } else {
        stop("Invalid estimator variant provided.")
    }
}

#' @title Estimate the N_pop size using a doubly robust method
#'
#' @description The doubly robust method estimates the total population size by
#' combining capture-recapture data with covariate information to improve
#' estimation accuracy.
#'
#' @details The doubly robust method leverages both the capture-recapture data
#' and covariate information to provide more accurate estimates of the total
#' population size. It uses cross-validation to estimate nuisance parameters
#' and combines these estimates to produce the final population size estimate.
#' It consists of 3 separate estimation techniques: DR, TMLE, and PI.
#'
#' @param data a data frame containing the observed capture histories and
#' covariates.
#' @param n_capture_cols the number of capture columns in the data.
#' @param estimator_variant an object specifying the doubly robust estimation
#' variant.
#' @return A list containing the estimation results and the total population
#' size.
#'
#' @importFrom drpop popsize
#' @export
estimate_doubly_robust <- function(
    data,
    n_capture_cols,
    estimator_variant = doubly_robust_estimator()
) {
    if (!inherits(estimator_variant, "doubly_robust_estimator")) {
        stop("Invalid estimator variant provided.")
    }

    n_folds <- estimator_variant$n_folds

    # Sufficient data is considered 50 observations per fold
    if (n_folds > 1 && n_folds > nrow(data) / 50) {
        n_folds <- pmax(floor(nrow(data) / 50), 1)
        warning(paste(
            "n_folds is reduced to",
            n_folds,
            "to have sufficient test data."
        ))
    }

    capture_cols <- paste0("capture_", seq_len(n_capture_cols))

    estimates <- data |>
        all_int_cols_to_numeric() |>
        extract_captured_data(capture_cols) |>
        drpop::popsize(
            K = n_capture_cols,
            funcname = estimator_variant$estimator_func,
            nfolds = n_folds,
            margin = estimator_variant$margin
        )
    return(estimates$result)
}

#' @title Estimate using Log Linear Methods
#'
#' @description Calculate an estimate for the total unobserved population by
#' fitting log linear models to the observed data.
#'
#' @details The classic CRC literature focuses on fitting log-linear models to
#' independent data and estimating the total population size based on the
#' observed capture histories. This function allows users to either estimate
#' the unobserved population by fitting one of two models (either `poisson` or
#' `negbin`) to the observed data and then selecting either via the minimized
#' AIC or a stepwise model selection procedure.
#'
#' @param data A data frame containing the observed capture histories.
#' @param n_capture_cols The number of capture columns in the data.
#' @param estimator_variant An object specifying the log linear estimator
#' variant to use.
#' @return A data frame containing the estimated total population size, the
#' lower and upper confidence intervals, and the AIC value if applicable.
#'
#' @examples
#' # Assuming `observed_data` is a data frame with capture columns
#' estimate_log_linear(observed_data, n_capture_cols = 3)
#'
#' @importFrom Rcapture
#' @export
estimate_log_linear <- function(
    data,
    n_capture_cols,
    estimator_variant = log_linear_estimator()
) {
    if (!inherits(estimator_variant, "log_linear_estimator")) {
        stop("Invalid estimator variant provided.")
    }

    if (estimator_variant$selection == "aic") {
        best_est <- aic_selection(data, n_capture_cols, estimator_variant) |>
            dplyr::mutate(
                estimate = estimate,
                cil = lower_ci,
                ciu = upper_ci,
                aic = AIC
            ) |>
            dplyr::select(estimate, cil, ciu, aic) |>
            dplyr::arrange(aic, ascending = TRUE) |>
            dplyr::filter(row_number() == 1)
        return(best_est)
    } else if (estimator_variant$selection == "stepwise") {
        return(stepwise_selection(data, n_capture_cols, estimator_variant))
    } else {
        stop("Invalid selection method provided.")
    }
}
