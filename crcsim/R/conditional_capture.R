################################################################################
# File: conditional_capture.R                                                  #
# Project: crcsim                                                              #
# Created Date: 2026-08-27                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################
#' Generate the conditional capture probability function.
#'
#' @param q_1 The nuisance function q_1.
#' @param q_2 The nuisance function q_2.
#' @param q_12 The nuisance function q_12.
#'
#' @return A function that computes the conditional capture probability given
#' covariate data.
#'
#' @references Das, Manjari and Kennedy, Edward H. and Jewell, Nicholas P.
#' (2024). "Doubly Robust Capture-Recapture Methods for Estimating Population
#' Size". Journal of the American Statistical Association, 119(546), 1309-1321.
#'
#' @export
conditional_capture <- function(q_1, q_2, q_12) {
    zero_probability <- !is.na(q_1) &
        !is.na(q_2) &
        (q_1 == 0 | q_2 == 0)
    if (any(zero_probability)) {
        warning("q_1 or q_2 is zero, which leads to division by zero.")
    }
    result <- q_12 / (q_1 * q_2)
    result[zero_probability] <- NA_real_
    return(result)
}
