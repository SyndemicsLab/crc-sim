################################################################################
# File: simulate.R                                                             #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-05                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Filter Contingency Table to Observed Rows. This returns only rows where all
#' capture variables are 0. Thus, the "unknown" rows are removed and the rows
#' with at least one observation are kept.
#'
#' @param data_table data.table: the contingency table to filter
#' @param capture list: strings of captures
#'
#' @keywords internal

filter_observed_rows <- function(data_table, capture) {
    data_table <- tibble::as_tibble(data_table)
    observed_rows <- data_table |>
        dplyr::mutate(
            tmp_capture_sum = rowSums(dplyr::across(dplyr::all_of(capture)))
        ) |>
        dplyr::filter(.data$tmp_capture_sum != 0) |>
        dplyr::select(-.data$tmp_capture_sum)

    return(observed_rows)
}

#' Apply Suppression to Frequency Column
#'
#' @keywords internal
#' @noRd

apply_suppression_filter <- function(data_table, suppress) {
    if (!isTRUE(suppress > 0)) {
        return(data_table)
    }

    return(data_table |> dplyr::filter(!(.data$N_ID %in% seq_len(suppress))))
}

#' Run CRC Model for Simulated Data
#'
#' @keywords internal
#' @noRd

run_simulation_model <- function(
    data_table,
    capture,
    method,
    formula_selection,
    opts_stepwise,
    opts_tmle
) {
    return(crc(
        data = data_table,
        freq_column = "N_ID",
        binary_variables = capture,
        method = method,
        formula_selection = formula_selection,
        opts_stepwise = opts_stepwise,
        opts_tmle = opts_tmle
    ))
}

#' Split Data by Group for Simulation
#'
#' @keywords internal
#' @noRd

split_simulation_groups <- function(data_table, group) {
    group_values <- sort(unique(data_table[[group]]))

    return(purrr::map(
        group_values,
        function(group_value) {
            group_table <- data_table |>
                dplyr::filter(.data[[group]] == group_value) |>
                dplyr::select(-dplyr::all_of(group))

            return(group_table)
        }
    ))
}

#' Attempts recovery of the ground truth - fundamentally a wrapper for crc
#'
#' @param data_table data.table from the \code{create.data} step
#' @param capture list: strings of captures
#' @param group stratification string to group by
#' @param suppress numeric: maximum value to suppress at
#' @param method string: selection for spatial capture recapture model - either
#' "poisson" or "negbin"
#' @param formula_selection string: selection for formula decision - either
#' "aic", "corr", "stepwise", or "tmle"
#' @param opts_stepwise list: list containing 'direction' of 'forward'
#' 'backward' or 'both', and 'threshold': p value threshold for stepwise
#' selection
#' @param opts_tmle list: TMLE options passed into \code{crc}
#'
#' @keywords internal

simulate <- function(
    data_table,
    group,
    suppress,
    method,
    formula_selection,
    opts_stepwise,
    opts_tmle,
    capture = c("APCD", "BSAS", "Casemix", "Death", "Matris", "PMP")
) {
    data_table <- filter_observed_rows(data_table, capture)
    data_table <- apply_suppression_filter(data_table, suppress)

    if (!missing(group)) {
        split_tables <- split_simulation_groups(data_table, group)
        out_list <- purrr::map(
            split_tables,
            run_simulation_model,
            capture = capture,
            method = method,
            formula_selection = formula_selection,
            opts_stepwise = opts_stepwise,
            opts_tmle = opts_tmle
        )

        return(out_list)
    }

    out <- run_simulation_model(
        data_table = data_table,
        capture = capture,
        method = method,
        formula_selection = formula_selection,
        opts_stepwise = opts_stepwise,
        opts_tmle = opts_tmle
    )

    return(out)
}
