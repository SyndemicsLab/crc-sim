################################################################################
# File: crossfit_fold.R                                                        #
# Project: crc-sim                                                             #
# Created Date: 2026-05-20                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-21                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Create cross-validation fold assignments for CRC estimation
#'
#' @param data A data frame containing the capture data. This should be the
#' raw' capture data, not a frequency table.
#' @param nfolds Integer scalar giving the number of cross-validation folds to
#' create. Must be a positive integer.
#' @return An integer vector of the same length as the number of rows in
#' \code{data}, where each element is an integer from 1 to \code{nfolds}
#' indicating the fold assignment for that row.
#'
#' @export
crossfit_fold <- function(data, nfolds) {
    if (!is.data.frame(data)) {
        stop("`data` must be a data frame.")
    }

    if (length(nfolds) != 1 || !is.numeric(nfolds) || nfolds < 1) {
        stop("`nfolds` must be a single positive number.")
    }

    nfolds <- as.integer(nfolds)
    n_rows <- nrow(data)

    permutset <- sample(1:n_rows, n_rows, replace = FALSE)
    idfold <- rep(1, n_rows)

    if (nfolds > 1) {
        test_indices <- lapply(seq_len(nfolds), function(folds) {
            sbset <- ((folds - 1) * ceiling(n_rows / nfolds) + 1):(folds *
                ceiling(n_rows / nfolds)) # nolint: indentation_linter
            sbset <- sbset[sbset <= n_rows]
            return(permutset[sbset])
        })

        idfold[unlist(test_indices, use.names = FALSE)] <- rep(
            seq_len(nfolds),
            lengths(test_indices)
        )
    }

    return(idfold)
}


build_list_pairs <- function(capture_columns) {
    if (length(capture_columns) < 2) {
        return(list())
    }
    combinations <- utils::combn(capture_columns, 2, simplify = FALSE)
    str_list <- unlist(lapply(combinations, paste, collapse = ","))
    return(list(combos = combinations, string_names = str_list))
}
