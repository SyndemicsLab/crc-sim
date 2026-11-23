################################################################################
# File: data_filtering.R                                                       #
# Project: crcsim                                                              #
# Created Date: 2026-10-21                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-21                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Apply Suppression to Frequency Table Data.
#' This function takes a frequency table dataset and applies suppression to the
#' frequency column.
#'
#' @param data data.frame: the frequency table dataset to apply suppression to
#' @param suppression_threshold int: the threshold at which to apply
#' suppression.
#' @param freq_col string: the name of the frequency column in the dataset.
#' Default is "N_ID".
#' @return a data.frame with suppression applied to the frequency column
#'
#' @importFrom dplyr mutate
#' @export
suppress_data <- function(data, suppression_threshold, freq_col = "N_ID") {
    if (!is.data.frame(data)) {
        stop("data must be a data.frame")
    }
    if (
        !is.numeric(suppression_threshold) || length(suppression_threshold) != 1
    ) {
        stop("suppression_threshold must be a numeric scalar")
    }
    if (!is.character(freq_col) || length(freq_col) != 1) {
        stop("freq_col must be a single column name")
    }
    if (!freq_col %in% colnames(data)) {
        stop("Frequency column not found in data. Cannot apply suppression.")
    }
    if (!is.numeric(data[[freq_col]])) {
        stop("Frequency column must be numeric")
    }
    data <- dplyr::filter(
        data,
        .data[[freq_col]] > suppression_threshold
    )
    return(data)
}

#' Get Captured Rows
#' @description Function to filter out rows where all capture columns are not 0.
#' This is used to identify all the captured rows in a dataset.
#'
#' @param data_table data.frame: Data frame containing either a frequency table
#' or individual records with capture columns.
#' @param capture character vector: names of capture columns
#' @returns tibble containing only rows where people are captured.
#'
#' @importFrom dplyr filter across
#' @importFrom tibble as_tibble
#' @export
extract_captured_data <- function(data_table, capture) {
    captured_rows <- tibble::as_tibble(data_table) |>
        dplyr::filter(rowSums(dplyr::across(dplyr::all_of(capture))) != 0)
    return(captured_rows)
}

#' Get Uncaptured Rows
#' @description Function to filter out rows where all capture columns are 0.
#' This is used to identify all the uncaptured rows in a dataset, which is
#' necessary for extracting a ground truth value for simulation evaluation.
#'
#' @param data_table data.frame: Data frame containing either a frequency table
#' or individual records with capture columns.
#' @param capture character vector: names of capture columns
#' @returns tibble containing only rows where people are uncaptured.
#'
#' @importFrom dplyr filter across
#' @importFrom tibble as_tibble
#' @export
extract_uncaptured_data <- function(data_table, capture) {
    uncaptured_rows <- tibble::as_tibble(data_table) |>
        dplyr::filter(rowSums(dplyr::across(dplyr::all_of(capture))) == 0)
    return(uncaptured_rows)
}
