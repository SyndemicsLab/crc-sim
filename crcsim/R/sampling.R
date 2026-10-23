################################################################################
# File: sampling.R                                                             #
# Project: crcsim                                                              #
# Created Date: 2026-10-21                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-21                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' @title Draw Binomial Samples
#'
#' @description A function to generate binomial samples along the specified
#' probability vector.
#'
#' @details It takes one required and an optional second parameter and returns
#' a vector of 0 or 1 integers. The vector is the same size as the number of
#' values in the `prob` vector.
#'
#' @param prob numeric vector: probabilities respective to each `n`.
#' @param prefix string: an optional prefix to append to each column name.
#' @return A named vector of 0 or 1 integers.
#'
#' @importFrom stats rbinom
#' @keywords internal
#' @noRd
binom_sample <- function(prob, prefix = "") {
    n <- length(prob)
    out <- as.integer(stats::rbinom(n = n, size = 1, prob = prob))
    names(out) <- paste0(prefix, 1:n)
    return(out)
}

#' @title Sample a Distribution data.frame
#'
#' @description This function takes in a distributional data.frame and samples
#' from it returning a list of sampled values.
#'
#' @details First the function checks that the data.frame is not empty and that
#' the appropriate distribution data.frame columns exist. If both succeed, it
#' draws a sample for each row in the data.frame and returns a list of sampled
#' values according to the distribution parameters.
#'
#' @param dist_df data.frame: A data.frame containing the distribution params
#' @param prefix string: A prefix to append to the sample names
#' @return A list containing the samples corresponding to the data.frame
#'
#' @keywords internal
#' @noRd
dist_df_sample <- function(dist_df, prefix = "") {
    n <- nrow(dist_df)
    if (n == 0) {
        return(list())
    }
    dist_df <- validate_dist_df(dist_df)
    samples <- purrr::pmap(dist_df, dist_sample)
    names(samples) <- paste0(prefix, seq_len(n))
    return(samples)
}

#' @title Sample a single distribution
#'
#' @description A helper function to sample a single distribution.
#'
#' @details Given the name of the distribution, the two governing parameters,
#' and the resulting data type, sample the value. It is generally used as a
#' helper function inside \code{dist_df_sample()}
#'
#' @param distribution string: The name, either `uniform` or `normal`
#' @param p1 numeric: First governing parameter
#' @param p2 numeric: Second governing parameter
#' @param dtype string: The data type to return
#' @return The sample for the distribution
#'
#' @importFrom stats runif rnorm
#' @keywords internal
#' @noRd
dist_sample <- function(distribution, p1, p2, dtype) {
    if (distribution == "uniform") {
        value <- stats::runif(n = 1, min = p1, max = p2)
    } else if (distribution == "normal") {
        value <- stats::rnorm(n = 1, mean = p1, sd = p2)
    } else {
        stop("Unsupported distribution in covariate specification")
    }

    if (dtype == "integer") {
        return(as.integer(round(value)))
    }

    out <- as.numeric(value)
    return(out)
}
