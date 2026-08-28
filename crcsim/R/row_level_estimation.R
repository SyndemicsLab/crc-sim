################################################################################
# File: row_level_estimation.R                                                 #
# Project: crcsim                                                              #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' CRC with Plugin Estimators
#' @description This function forms the basis for the CRC estimation process on
#' row level data. This usually happens when there are covariates for each
#' observation that we can use during the estimation process. The function
#' takes in a data frame and an \code{EstimatorOptions} object specifying the
#' plugin estimator to use and its parameters. The function first checks that
#' the provided options are valid and that the data is not a frequency table.
#' It then identifies the capture indicator columns, either from the options or
#' by inferring them from the data. Next, it generates cross-validation folds
#' and computes capture probabilities using the specified nuisance function and
#' plugin estimator. Finally, it calculates and returns the TMLE estimate of
#' population size along with confidence intervals. Supported plugin estimators
#' include "logit", "ranger", and "gam".
#'
#' @param data a data frame containing the capture data. All columns must be
#' numeric type as long as we use the \code{drpop} package.
#' @param opts a \code{EstimatorOptions} object specifying the plugin method and
#'  parameters for the TMLE estimation.
#' @returns a list containing the estimate of population size and
#' confidence intervals.
#'
#' @importFrom drpop popsize
#' @export
row_level_estimation <- function(data, opts) {
    if (!inherits(opts, "EstimatorOptions")) {
        stop(paste(
            "Invalid EstimatorOptions object provided.",
            "Cannot run CRC with plugin estimator."
        ))
    }
    if (is_frequency_table(data, "N_ID")) {
        stop(paste(
            "Data must not be a frequency table. ",
            "Please provide raw capture data for plugin estimation."
        ))
    }

    if (is.null(opts[["capture_columns"]])) {
        capture_columns <- names(
            data[, vapply(data, is_column_binary, logical(1L))]
        )
    } else {
        capture_columns <- opts[["capture_columns"]]
    }
    sim_data <- data |>
        all_int_cols_to_numeric() |>
        extract_captured_data(capture_columns)

    n <- nrow(sim_data)

    estimates <- estimate_capture_prob(
        sim_data,
        length(capture_columns),
        method = opts[["model"]],
        func = opts[["nusiance_function"]],
        nfolds = opts[["nfolds"]],
        margin = opts[["threshold"]]
    )

    return(estimates)
}
