################################################################################
# File: row_level_estimation.R                                                 #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-06-16                                                    #
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

    # This is the basic outline of what we want to do rather than utilize drpop.
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
    n_folds <- opts[["nfolds"]]
    funcname <- opts[["model"]]
    margin <- opts[["threshold"]]
    qhat <- drpop::popsize(
        data = sim_data,
        funcname = funcname,
        nfolds = n_folds,
        margin = margin
    )

    psin_estimates <- drpop::popsize(
        data = sim_data,
        getnuis = qhat[["nuis"]],
        idfold = qhat[["idfold"]]
    )
    return(psin_estimates)

    # folds <- crossfit_fold(data, opts[["nfolds"]])

    # cap_probs <- lapply(seq_len(opts[["nfolds"]]), function(fold) {
    #     test_idx <- folds == fold
    #     train_idx <- !test_idx
    #     capture_prob <- qhat_generation(
    #         train = data[train_idx, , drop = FALSE],
    #         test = data[test_idx, , drop = FALSE],
    #         capture_names = capture_columns,
    #         nuisance_function = opts[["model"]]
    #     ) |>
    #         estimate_capture_probability(opts[["estimator"]])
    #     return(capture_prob)
    # })

    # return(list(
    #     estimate = get_n_estimate(capture_prob),
    #     lower_ci = get_lower_ci(capture_prob),
    #     upper_ci = get_upper_ci(capture_prob),
    #     estimator = opts$estimator
    # ))
}
