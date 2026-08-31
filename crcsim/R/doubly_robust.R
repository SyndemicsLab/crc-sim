################################################################################
# File: doubly_robust.R                                                        #
# Project: crcsim                                                              #
# Created Date: 2026-08-27                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Compute the doubly robust estimator for capture-recapture data.
#'
#' @param q_1 The nuisance function q_1.
#' @param q_2 The nuisance function q_2.
#' @param q_12 The nuisance function q_12.
#' @param phi The auxiliary function phi.
#'
#' @return The doubly robust estimate of the capture-recapture parameter.
#'
#' @export
doubly_robust_estimator <- function(q_1, q_2, q_12, phi) {
    plugin_phi <- plugin_estimator(q_1, q_2, q_12)
    phi_over_q <- mean(phi, na.rm = TRUE)
    psi_inverse_hat <- max(plugin_phi + phi_over_q, 1)
    return(psi_inverse_hat)
}
