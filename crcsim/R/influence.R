################################################################################
# File: influence.R                                                            #
# Project: crcsim                                                              #
# Created Date: 2026-08-27                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

# Build the influence function. This should require the conditional capture probability function, the nuisance functions (q1, q2, and q12), datasets corresponding to Y_1 and Y_2, and the conditional capture probability

#' Generate the influence function for the given data and nuisance functions.
#'
#' @param covariate_data The covariate data for the observations.
#' @param list_1 The dataset corresponding to Y_1.
#' @param list_2 The dataset corresponding to Y_2.
#' @param q_1 The nuisance function q_1.
#' @param q_2 The nuisance function q_2.
#' @param q_12 The nuisance function q_12.
#' @param psi The nuisance function psi.
#'
#' @return The influence function evaluated for the given data and nuisance
#' functions.
#'
#' @references Das, Manjari and Kennedy, Edward H. and Jewell, Nicholas P.
#' (2024). "Doubly Robust Capture-Recapture Methods for Estimating Population
#' Size". Journal of the American Statistical Association, 119(546), 1309-1321.
influence <- function(
    list_1,
    list_2,
    q_1,
    q_2,
    q_12
) {
    gamma_inv <- 1 / conditional_capture(q_1, q_2, q_12)
    psi_inv <- plugin_estimator(q_1, q_2, q_12)

    bias_correction <- (list_1 / q_1) +
        (list_2 / q_2) -
        (list_1 * list_2 / q_12)

    phi <- gamma_inv * bias_correction - psi_inv
    return(phi)
}
