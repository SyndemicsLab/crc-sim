################################################################################
# File: ranger.R                                                               #
# Project: crc-sim                                                             #
# Created Date: 2026-05-14                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' TMLE using ranger
#' This function performs TMLE estimation of population size using ranger
#' as the nuisance function estimator. It is currently a wrapper for the drpop
#' package.
#'
#' @param data a data frame containing the capture data. All columns must be
#' numeric type as long as we use the \code{drpop} package.
#' @param nfolds int: number of folds to use for cross-fitting. Default is 2.
#' @param margin numeric: margin parameter for the TMLE estimation. Default is
#' 0.005.
#' @returns a list containing the TMLE estimate of population size and
#' confidence intervals.
#'
#' @importFrom drpop popsize
#' @export
ranger_estimate <- function(data, nfolds = 2, margin = 0.005) {
    qhat <- data |>
        all_int_cols_to_numeric() |>
        popsize(
            funcname = "ranger",
            nfolds = nfolds,
            margin = margin
        )
    psin_estimate <- popsize(
        data = data,
        getnuis = qhat$nuis,
        idfold = qhat$idfold
    )
    return(psin_estimate)
}
