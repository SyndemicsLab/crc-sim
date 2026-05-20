################################################################################
# File: extract.R                                                              #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-19                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Filter Captured Rows
#' @description Internal function to filter rows where all capture columns
#' are 0.
#'
#' @param data_table data.frame: contingency table containing capture columns
#' @param capture character vector: names of capture columns
#' @returns tibble containing only uncaptured rows
#' @keywords internal
#' @noRd
filter_captured_rows <- function(data_table, capture) {
    captured_rows <- tibble::as_tibble(data_table) |>
        dplyr::filter(rowSums(across(capture)) != 0)
    return(captured_rows)
}

#' Filter Uncaptured Rows
#' @description Internal function to filter rows where all capture columns
#' are 0.
#'
#' @param data_table data.frame: contingency table containing capture columns
#' @param capture character vector: names of capture columns
#' @returns tibble containing only uncaptured rows
#' @keywords internal
#' @noRd
filter_uncaptured_rows <- function(data_table, capture) {
    uncaptured_rows <- tibble::as_tibble(data_table) |>
        dplyr::filter(rowSums(across(capture)) == 0)
    return(uncaptured_rows)
}

#' Build Uncaptured Count Values
#'
#' @param uncaptured_rows data.frame: rows where all capture columns are 0
#' @returns list of uncaptured frequency counts
#' @keywords internal
#' @noRd
build_uncaptured_count_values <- function(uncaptured_rows) {
    return(as.list(uncaptured_rows[["N_ID"]]))
}

#' Build Uncaptured Count Names
#'
#' @param uncaptured_rows data.frame: rows where all capture columns are 0
#' @param group character: strata column name used for naming output
#' @param group_missing logical: whether the group argument was omitted
#' @returns character vector or list used as output names
#' @keywords internal
#' @noRd
build_uncaptured_count_names <- function(
    uncaptured_rows,
    group,
    group_missing
) {
    if (group_missing) {
        return(rep("base", nrow(uncaptured_rows)))
    }

    return(as.list(uncaptured_rows[[group]]))
}

#' Extract Uncaptured Ground Truth
#'
#' Extracts the uncaptured count from a contingency table where uncaptured
#' observations are defined as rows with all capture columns equal to 0.
#'
#' @param data_table data.frame: contingency table from \code{create_data}
#' @param group character: strata column used to name output values
#' @param capture character vector: names of capture columns
#' @returns named list of uncaptured frequency counts
#' @export
extract_ground_truth <- function(
    data_table,
    group,
    capture = c("APCD", "BSAS", "Casemix", "Death", "Matris", "PMP")
) {
    group_missing <- missing(group)
    uncaptured_rows <- filter_uncaptured_rows(data_table, capture)
    out <- build_uncaptured_count_values(uncaptured_rows)
    out_names <- build_uncaptured_count_names(
        uncaptured_rows,
        group,
        group_missing
    )

    names(out) <- out_names

    return(out)
}
