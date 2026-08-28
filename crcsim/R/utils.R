################################################################################
# File: utils.R                                                                #
# Project: crcsim                                                              #
# Created Date: 2026-08-28                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-28                                                    #
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
    list1_vec <- 1:(n_lists - 1)
    list2_vec <- 1:n_lists
    result <- tidyr::expand_grid(j1 = list1_vec, k1 = list2_vec) |>
        dplyr::filter(k1 > j1) |>
        dplyr::mutate(pair = paste0(j1, ",", k1), .keep = "none") |>
        dplyr::pull(pair)
    return(result)
}

#' Test the overlap between two capture lists.
#'
#' @param l1 The first capture list.
#' @param l2 The second capture list.
#' @param margin The minimum acceptable overlap between the lists.
#' @return NULL. A warning is issued if the overlap is less than the margin.
#'
#' @export
test_list_overlap <- function(l1, l2, margin) {
    if (mean(l1 * l2) < margin) {
        warning(paste0(
            "Warning: Overlap between the lists ",
            j,
            " and ",
            k,
            " is less than ",
            margin,
            ".\n"
        ))
    }
    return(NULL)
}

#' Validate the capture columns contain only 0 or 1 values
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
        data[, start:end],
        2,
        function(x) all(x %in% c(0, 1))
    )
    n_lists <- end - start + 1
    if (sum(zero_one_check) != n_lists) {
        return(FALSE)
    }
    return(TRUE)
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
