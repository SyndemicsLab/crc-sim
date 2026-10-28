################################################################################
# File: utils.R                                                                #
# Project: crcsim                                                              #
# Created Date: 2026-08-28                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-27                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Build all unique pairs of capture lists.
#'
#' @param n_lists The number of capture lists.
#' @return A character vector containing all unique pairs of capture lists in
#' the format "j,k".
#'
#' @importFrom tidyr expand_grid
#' @importFrom dplyr filter mutate pull
#' @export
build_listpairs <- function(n_lists) {
    if (n_lists < 2) {
        return(character())
    }

    list1_vec <- 1:(n_lists - 1)
    list2_vec <- 1:n_lists
    result <- tidyr::expand_grid(j1 = list1_vec, k1 = list2_vec) |>
        dplyr::filter(k1 > j1) |>
        dplyr::mutate(pair = paste0(j1, ",", k1), .keep = "none") |>
        dplyr::pull(pair)
    return(result)
}

#' Validate if a value is a numeric scalar.
#'
#' @param x The value to check.
#' @return TRUE if the value is a numeric scalar, FALSE otherwise.
#'
#' @export
is_valid_number <- function(x) {
    return(is.numeric(x) && length(x) == 1 && !is.na(x))
}

#' Test the overlap between two capture lists.
#'
#' @param l1 The first capture list.
#' @param l2 The second capture list.
#' @param margin The minimum acceptable overlap between the lists.
#' @return NULL. A warning is issued if the overlap is less than the margin.
#'
#' @export
test_list_overlap <- function(l1, l2, margin, l1_name = NULL, l2_name = NULL) {
    overlap <- mean(l1 * l2)
    if (overlap < margin) {
        error_prefix <- paste0(
            "Overlap between ",
            ifelse(is.null(l1_name), "the first list", l1_name)
        )
        error_prefix <- paste0(
            error_prefix,
            " and ",
            ifelse(is.null(l2_name), "the second list", l2_name),
            " is "
        )
        warning(paste0(
            error_prefix,
            signif(overlap, 3),
            " less than the threshold ",
            signif(margin, 3),
            "."
        ))
    }
    return(NULL)
}

#' Create a zero matrix with specified dimensions and optional row and column
#' names.
#'
#' @param nrow The number of rows in the matrix.
#' @param cols The column names for the matrix.
#' @param rownames The row names for the matrix (optional).
#' @return A matrix filled with zeros with the specified dimensions and names.
#'
#' @export
zero_matrix <- function(nrow, cols, rownames = NULL) {
    zero_mat <- matrix(
        0,
        nrow = nrow,
        ncol = length(cols),
        dimnames = list(rownames, cols)
    )
    return(zero_mat)
}

#' Get the names of the capture columns in the dataset. Note: this assumes all capture columns start with the prefix `capture_`.
#'
#' @param data The dataset to inspect.
#' @return A character vector of capture column names.
#'
#' @export
get_capture_cols <- function(data) {
    return(grep("^capture_", names(data), value = TRUE))
}

#' Get the observed data based on capture columns.
#'
#' @param data The dataset to filter.
#' @param capture_columns The capture columns to consider.
#' @return A subset of the data where at least one capture column is positive.
#'
#' @export
get_observed_data <- function(data, capture_columns) {
    return(data[rowSums(data[capture_columns]) > 0, , ])
}

#' Get the unobserved data based on capture columns. This will return all rows
#' where the `capture_columns` are all zero.
#'
#' @param data The dataset to filter.
#' @param capture_columns The capture columns to consider.
#' @return A subset of the data where all capture columns are zero.
#'
#' @export
get_unobserved_data <- function(data, capture_columns) {
    return(data[rowSums(data[capture_columns]) == 0, , ])
}
