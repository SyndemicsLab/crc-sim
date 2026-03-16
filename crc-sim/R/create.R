################################################################################
# File: create.R                                                               #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-02-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Data generation tool
#'
#' @param n_individuals int: number of participants to simulate
#' @param p_captures list: list of capture probabilities to emulate
#' @param p_strata list: list of extra 'stratification' probabilities to emulate
#'
#' @import data.table
#' @returns a data.table
#'
#' @export
create_data <- function(n_individuals, p_captures, p_strata) {
    data_table <- lapply(1:n_individuals, function(x) {
        return(
            data.table(
                t(create_capture(length(p_captures), p_captures)),
                strata = create_strata(length(p_strata), p_strata)
            )
        )
    }) |>
        data.table::rbindlist()
    return(data_table[, list(N_ID = .N), by = c(names(data_table))])
}

#' A function to sample from 1:n along probability \code{prob}
#'
#' @param n int: number of 'strata' to simulate
#' @param prob list: probabilities of n - expects summation to 1
#'
#' @keywords internal
create_strata <- function(n, prob) {
    if (length(prob) != n) {
        stop("Probability length differs from n")
    }
    if (sum(prob) != 1) {
        warning("Probability does not sum to 1")
    }

    return(sample(1:n, 1, prob = prob))
}

#' A function to generate binomial flips along the probability specified
#'
#' @param n int: number of 'captures' to simulate
#' @param prob list: probabilities respective to each \code{n}
#'
#' @importFrom stats rbinom
#' @keywords internal
create_capture <- function(n, prob) {
    if (length(prob) != n) {
        stop("Probability length differs from n")
    }

    out <- vector(length = n)
    for (i in seq_along(1:n)) {
        out[i] <- rbinom(1, 1, prob[[i]])
    }
    names(out) <- paste0("capture_", 1:n)

    return(out)
}
