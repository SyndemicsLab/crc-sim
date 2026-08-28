################################################################################
# File: crossfold.R                                                            #
# Project: crcsim                                                              #
# Created Date: 2026-08-28                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Build cross-validation folds and run estimation for each fold.
#'
#' @param fold The current fold number.
#' @param data The dataset to be split into folds.
#' @param nfolds The total number of folds.
#' @param permutation A permutation of row indices for creating folds.
#' @return A list containing the training set, test set, and the subset of
#' indices for the current fold.
#'
#' @export
build_folds <- function(fold, data, nfolds, permutation) {
    if (nfolds == 1) {
        sbset <- seq_len(nrow(data))
        train <- data
        test <- data
    } else {
        fold_size <- ceiling(nrow(data) / nfolds)
        sbset <- ((fold - 1) * fold_size + 1):(fold * fold_size)
        sbset <- sbset[sbset <= nrow(data)]
        train <- data[permutation[-sbset], ]
        test <- data[permutation[sbset], ]
    }

    return(list(
        train = train,
        test = test,
        sbset = sbset
    ))
}

#' Run a single fold of the cross-validation procedure.
#'
#' @param fold_split A list containing the training set, test set, and the
#' subset of indices for the current fold.
#' @param method The estimation method to be used.
#' @param j The index of the first list.
#' @param k The index of the second list.
#' @param margin The margin parameter for estimation.
#' @param n_lists The total number of capture lists.
#' @param n The total number of observations.
#' @param nuisance_estimation_func The function used for nuisance parameter
#' estimation.
#' @param ... Additional arguments passed to the nuisance estimation function.
#' @return A list containing the capture probability and its variance for the
#' current fold.
#'
#' @export
run_fold <- function(
    fold_split,
    method,
    j,
    k,
    margin,
    n_lists,
    n,
    nuisance_estimation_func
) {
    train <- fold_split[["train"]]
    test <- fold_split[["test"]]
    test_list_overlap(train[, j], train[, k], margin)

    nuisance_functions <- try(
        nuisance_estimation(func, train, test, n_lists, j, k, margin),
        silent = TRUE
    )

    if (inherits(nuisance_functions, "try-error")) {
        return(NULL)
    }

    # Returns a list of capture_probability and variance
    return(estimate_psi(
        method,
        nuisance_functions,
        y_j = test[, j],
        y_k = test[, k],
        margin = margin,
        n_lists = n_lists,
        n = n
    ))
}
