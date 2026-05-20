################################################################################
# File: data_validation.R                                                      #
# Project: crc-sim                                                             #
# Created Date: 2026-05-14                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' A function used to check if a column in a data frame is binary
#' (i.e., contains only 0s and 1s).
#'
#' @param column a vector representing a column in a data frame
#' @returns a boolean indicating whether the column is binary
#'
#' @keywords internal
#' @export
is_column_binary <- function(column) {
    unique_values <- unique(column)
    return(length(unique_values) == 2 && all(unique_values %in% c(0, 1)))
}

#' A function used to check if the provided data is a frequency table. This is
#' necessary as the CRC estimation functions require a frequency table as input.
#'
#' @param data a data frame
#' @param frequency_column a string specifying the name of the frequency column
#' in the data frame
#' @returns a boolean indicating whether the data frame is a frequency table
#' (i.e., contains the specified frequency column and that column is numeric)
#'
#' @keywords internal
#' @export
is_frequency_table <- function(data, frequency_column) {
    if (
        frequency_column %in%
            names(data) &&
            is.numeric(data[[frequency_column]])
    ) {
        return(TRUE)
    }
    return(FALSE)
}

#' A function used to convert all integer columns to numeric. This is necessary
#' as the drpop package requires all data columns to be numerics.
#' This is a temporary workaround until we can implement our own TMLE estimation
#' function that can handle integer columns.
#'
#' @param df a data frame
#' @returns a data frame with all integer columns converted to numeric
#'
#' @importFrom dplyr mutate across where
#' @keywords internal
#' @export
all_int_cols_to_numeric <- function(df) {
    return(mutate(df, across(where(is.integer), as.numeric)))
}
