################################################################################
# File: validators.R                                                           #
# Project: crcsim                                                              #
# Created Date: 2026-10-21                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Validate the population size in the simulation.
#'
#' @param n_pop Number of individuals to simulate.
#' @return Invisibly returns NULL if the value is valid; otherwise, throws an
#'   error.
#'
#' @keywords internal
#' @noRd
validate_n_pop <- function(n_pop) {
    if (!is.numeric(n_pop) || length(n_pop) != 1) {
        stop("n_pop must be a numeric scalar")
    }
    if (n_pop <= 0 || n_pop %% 1 != 0) {
        stop("n_pop must be a positive integer")
    }
    return(invisible(NULL))
}

#' Validate the number of capture datasets.
#'
#' @param n_capture_cols Number of capture datasets and capture columns.
#' @return Invisibly returns NULL if the value is valid; otherwise, throws an
#'   error.
#'
#' @keywords internal
#' @noRd
validate_n_capture_cols <- function(n_capture_cols) {
    if (
        !is_valid_number(n_capture_cols) ||
            n_capture_cols < 1 ||
            n_capture_cols %% 1 != 0
    ) {
        stop("n_capture_cols must be a positive integer")
    }
    return(invisible(NULL))
}

#' @title Validate a Probability Vector
#'
#' @description Function to validate that a given vector contains valid
#' probabilities between 0 and 1.
#'
#' @param probs A vector of probabilities to validate.
#' @return Invisibly returns NULL if the vector is valid; otherwise, throws an
#' error.
#'
#' @keywords internal
validate_prob_vector <- function(probs) {
    if (!is.numeric(probs) || length(probs) == 0 || anyNA(probs)) {
        stop("probs must be a non-empty numeric vector without NA")
    }
    if (any(probs < 0 | probs > 1)) {
        stop("probs must be between 0 and 1")
    }
    return(invisible(NULL))
}

#' @title Validate Inputs for create_data.
#'
#' @description This checks that all inputs to the \code{create_data} function
#' are of the correct type and value range. If any checks fail, an error is
#' thrown.
#'
#' @details Exclusively used within \code{create_data} and not for public
#' consumption. This function verifies that `n_pop` and `n_capture_cols` are
#' valid numeric scalar positive integers, `cov_dists` is either NULL or a
#' valid \code{data.frame}, and that `return_frequency_table` is a logical
#' scalar. As this is a validation function with a lot of if checks, we ignore
#' the cyclocomp_linter.
#'
#' @param n_pop The number of people to simulate in the dataset.
#' @param n_capture_cols The number of capture columns.
#' @param cov_dists The possible range of covariate data.
#' @param freq_table A flag to indicate if a frequency table should be returned.
#' @return NULL on success or stops program with an error
#'
#' @keywords internal
#' @noRd
validate_create_data_params <- function(
    n_pop,
    n_capture_cols,
    cov_dists,
    freq_table
) {
    validate_n_pop(n_pop)
    validate_n_capture_cols(n_capture_cols)
    if (!is.null(cov_dists) && !is.data.frame(cov_dists)) {
        stop("cov_dists must be a data.frame or NULL")
    }
    if (!is.logical(freq_table) || length(freq_table) != 1) {
        stop("freq_table must be a logical scalar")
    }

    return(invisible(NULL))
}

#' @title Validate Binary Columns
#'
#' @description Validate the provided columns contain only 0 or 1 values.
#'
#' @param data The dataset containing capture columns.
#' @param end The number of capture columns to validate.
#' @param start The starting column index for the capture columns.
#' @return NULL if all capture columns contain only 0 or 1 values; otherwise,
#' an error is thrown.
#'
#' @keywords internal
validate_binary_cols <- function(data, end, start = 1) {
    zero_one_check <- apply(
        data[, start:end, drop = FALSE],
        2,
        function(x) all(x %in% c(0, 1))
    )
    n_lists <- end - start + 1
    if (sum(zero_one_check) != n_lists) {
        return(FALSE)
    }
    return(TRUE)
}

#' @title Validate Capture Probabilities.
#'
#' @description This first checks that the capture probabilities and handles
#' them accordingly returning a valid capture probability set.
#'
#' @details In handling the capture probabilities, it checks that they're
#' specified and if not returns the default \code{rep(0.5, n_capture_cols)}. If
#' probabilities are specified, this checks that they are numeric, of the
#' correct length, and between 0 and 1. If any checks fail, an error is thrown.
#' This is an internal function used within \code{create_data}.
#'
#' @param capture_probs The probabilities to verify against the capture cols.
#' @param n_capture_cols The capture columns to verify against the probs.
#' @return A valid capture probability vector or stops execution
#'
#' @keywords internal
#' @noRd
validate_capture_probs <- function(capture_probs, n_capture_cols) {
    if (is.null(capture_probs)) {
        return(rep(0.5, n_capture_cols))
    }

    if (!is.numeric(capture_probs)) {
        stop("capture_probs must be a numeric vector")
    }
    if (length(capture_probs) != n_capture_cols) {
        stop("Length of capture_probs must equal n_capture_cols")
    }
    if (any(capture_probs < 0 | capture_probs > 1)) {
        stop("All capture_probs values must be between 0 and 1")
    }

    return(as.numeric(capture_probs))
}

#' @title Validate Distribution `data.frame`
#'
#' @description This validates a distribution data.frame has correct columns.
#'
#' @details Distribution `data.frame` objects must have the following columns:
#' - `distribution`: A string of either `uniform` or `normal`
#' - `p1`: The first governing parameter of a distribution
#' - `p2`: The second governing parameter of a distribution
#' - `dtype`: The datatype of the parameter
#' And each row corresponds to a different distribution to sample from.
#'
#' @param dist_df data.frame: A data.frame object containing distribution data.
#' @return The validated data.frame or stops execution
#'
#' @keywords internal
#' @export
validate_dist_df <- function(dist_df) {
    # verify correct columns
    required_columns <- c("distribution", "p1", "p2", "dtype")
    if (!all(required_columns %in% names(dist_df))) {
        stop(
            paste(
                "A distribution dataframe must contain columns:",
                "distribution, p1, p2, dtype"
            )
        )
    }
    dist_df <- dist_df[, required_columns, drop = FALSE]
    return(dist_df)
}

#' @title Validate Covariate Specifications.
#'
#' @description This validates and standardizes \code{cov_dists}.
#'
#' @details In validating `cov_dists` we expect the following columns:
#' - `distribution`: The name of the distribution to draw from.
#' - `p1`: The first governing parameter of the distribution.
#' - `p2`: The second governing parameter of the distribution.
#' - `dtype`: The data type of the covariate.
#'
#' Each row in the `cov_dists` is treated such that it defines one generated
#' covariate. As this is a validation function with a lot of if checks, we
#' ignore cyclocomp_linter.
#'
#' @param cov_dists A dataframe describing the covariate information
#' @return A valid cov_dists dataframe or stops execution
#'
#' @keywords internal
#' @noRd
# nolint start: cyclocomp_linter
validate_cov_dists <- function(cov_dists) {
    if (is.null(cov_dists) || nrow(cov_dists) == 0) {
        return(data.frame(
            distribution = character(0),
            p1 = numeric(0),
            p2 = numeric(0),
            dtype = character(0)
        ))
    }

    cov_dists <- validate_dist_df(cov_dists)

    # string processing to ensure we know what we're looking at
    cov_dists[["distribution"]] <- tolower(as.character(
        cov_dists[["distribution"]]
    ))
    cov_dists[["dtype"]] <- tolower(as.character(cov_dists[["dtype"]]))

    # validate that distributions are either `uniform` or `normal`
    invalid_distributions <- setdiff(
        unique(cov_dists[["distribution"]]),
        c("uniform", "normal")
    )
    if (length(invalid_distributions) > 0) {
        stop("distribution must be one of: uniform, normal")
    }

    # validate that dtypes are either `numeric` or `integer`
    invalid_dtypes <- setdiff(
        unique(cov_dists[["dtype"]]),
        c("numeric", "integer")
    )
    if (length(invalid_dtypes) > 0) {
        stop("dtype must be one of: numeric, integer")
    }

    # validate distribution parameters are numbers
    if (!is.numeric(cov_dists[["p1"]]) || !is.numeric(cov_dists[["p2"]])) {
        stop("p1 and p2 must be numeric")
    }
    if (any(is.na(cov_dists[["p1"]]) | is.na(cov_dists[["p2"]]))) {
        stop("p1 and p2 cannot contain NA values")
    }

    # handle invalid uniform parameters
    uniform_rows <- cov_dists[["distribution"]] == "uniform"
    if (
        any(cov_dists[["p1"]][uniform_rows] > cov_dists[["p2"]][uniform_rows])
    ) {
        stop("For uniform, p1 (low) must be <= p2 (high)")
    }

    # handle invalid normal parameters
    normal_rows <- cov_dists[["distribution"]] == "normal"
    if (any(cov_dists[["p2"]][normal_rows] <= 0)) {
        stop("For normal, p2 (standard deviation) must be > 0")
    }

    return(cov_dists)
}
# nolint end: cyclocomp_linter

#' @title Validate and Set a Seed
#'
#' @description Verify the provided seed is either NULL or an integer scalar
#' and set the value for R random number generation.
#'
#' @details If the seed is NULL, do not set a seed. If the seed is an integer
#' scalar, we use that as a seed for the RNG engine.
#'
#' @param seed An integer scalar or NULL.
#' @return Invisibly returns NULL if the seed is valid; otherwise, throws an
#' error.
#'
#' @keywords internal
#' @noRd
# nolint start: cyclocomp_linter
validate_and_set_seed <- function(seed) {
    if (is.null(seed)) {
        return(invisible(NULL))
    }
    if (
        !is.numeric(seed) ||
            length(seed) != 1 ||
            is.na(seed) ||
            seed %% 1 != 0
    ) {
        stop("seed must be NULL or an integer scalar")
    }
    set.seed(seed)
    return(invisible(NULL))
}
# nolint end: cyclocomp_linter

#' @title Validate Covariate Data is Supplied
#'
#' @description Helper function to validate that covariate data, alpha, and
#' beta are correctly supplied and consistent with the number of captures. This
#' is used prior to attempting to create covariate based capture probabilities.
#'
#' @param n_capture_cols The number of capture columns expected in the
#' covariate data.
#' @param cov_dists The covariate distributions.
#' @param alpha The alpha parameter for the covariate model.
#' @param beta The beta parameter for the covariate model.
#' @return Invisibly returns NULL if the covariate data is valid; otherwise,
#' throws an error.
#'
#' @keywords internal
#' @noRd
# nolint start: cyclocomp_linter
validate_covariate_scenario <- function(
    n_capture_cols,
    cov_dists,
    alpha,
    beta
) {
    if (is.null(alpha) || is.null(beta)) {
        stop("alpha and beta must be supplied together")
    }
    if (is.null(cov_dists)) {
        stop("cov_dists is required with alpha and beta")
    }
    if (!is.numeric(alpha) || length(alpha) == 0 || anyNA(alpha)) {
        stop("alpha must be a non-empty numeric vector without NA")
    }
    if (!is.numeric(beta) || !is.matrix(beta) || anyNA(beta)) {
        stop("beta must be a numeric matrix without NA")
    }
    if (nrow(beta) != length(alpha)) {
        stop("beta must have one row for each alpha value")
    }
    if (ncol(beta) != nrow(cov_dists)) {
        stop("beta must have one column for each covariate")
    }
    if (!is_valid_number(n_capture_cols) || n_captures != length(alpha)) {
        stop("n_capture_cols must be equal to the length of alpha")
    }
    return(invisible(NULL))
}
# nolint end: cyclocomp_linter
