################################################################################
# File: create.R                                                               #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-05                                                    #
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
#' @returns a data.table
#'
#' @export
create_data <- function(n_individuals, p_captures, p_strata) {
    records <- build_individual_records(n_individuals, p_captures, p_strata)
    aggregated <- aggregate_contingency_table(records)

    return(data.table::as.data.table(aggregated))
}

#' Build Individual Simulation Records
#'
#' @keywords internal
#' @noRd
build_individual_records <- function(n_individuals, p_captures, p_strata) {
    records <- purrr::map_dfr(
        seq_len(n_individuals),
        create_individual_record,
        p_captures = p_captures,
        p_strata = p_strata
    )

    return(records)
}

#' Create One Individual Simulation Record
#'
#' @keywords internal
#' @noRd
create_individual_record <- function(index, p_captures, p_strata) {
    capture_values <- create_capture(length(p_captures), p_captures)
    strata_value <- create_strata(length(p_strata), p_strata)

    return(tibble::as_tibble_row(c(capture_values, strata = strata_value)))
}

#' Aggregate Records into a Contingency Table
#'
#' @keywords internal
#' @noRd
aggregate_contingency_table <- function(records) {
    output <- records |>
        dplyr::count(dplyr::across(dplyr::everything()), name = "N_ID")

    return(output)
}

#' A function to sample from 1:n along probability \code{prob}
#'
#' @param n int: number of 'strata' to simulate
#' @param prob list: probabilities of n - expects summation to 1
#'
#' @keywords internal
#' @noRd
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
#' @noRd
create_capture <- function(n, prob) {
    if (length(prob) != n) {
        stop("Probability length differs from n")
    }

    out <- stats::rbinom(n = n, size = 1, prob = unlist(prob))
    names(out) <- paste0("capture_", 1:n)

    return(out)
}
