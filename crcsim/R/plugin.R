################################################################################
# File: plugin.R                                                               #
# Project: crcsim                                                              #
# Created Date: 2026-08-27                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-27                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Compute the plugin estimator for the capture probability.
#'
#' @param q_1 The nuisance function q_1.
#' @param q_2 The nuisance function q_2.
#' @param q_12 The nuisance function q_12.
#'
#' @return The plugin estimate of the capture probability.
#'
#' @references Das, Manjari and Kennedy, Edward H. and Jewell, Nicholas P.
#' (2024). "Doubly Robust Capture-Recapture Methods for Estimating Population
#' Size". Journal of the American Statistical Association, 119(546), 1309-1321.
#'
#' @export
plugin_estimator <- function(q_1, q_2, q_12) {
    gamma_hat <- conditional_capture(q_1, q_2, q_12)
    gamma_inv_hat <- 1 / gamma_hat
    psi_inverse_hat <- mean(gamma_inv_hat)
    return(psi_inverse_hat)
}
