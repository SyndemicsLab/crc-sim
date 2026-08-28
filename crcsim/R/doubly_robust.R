################################################################################
# File: doubly_robust.R                                                        #
# Project: crcsim                                                              #
# Created Date: 2026-08-27                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-27                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

doubly_robust_estimator <- function(q_1, q_2, q_12, phi) {
    plugin_phi <- plugin_estimator(q_1, q_2, q_12)
    phi_over_q <- mean(phi, na.rm = TRUE)
    psi_inverse_hat <- max(plugin_phi + phi_over_q, 1)
    return(psi_inverse_hat)
}
