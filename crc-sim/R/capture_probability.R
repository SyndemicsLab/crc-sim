################################################################################
# File: capture_probability.R                                                  #
# Project: crc-sim                                                             #
# Created Date: 2026-05-18                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-20                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

conditional_capture_prob <- function(q_j, q_k, q_jk) {
    if (q_j == 0 || q_k == 0) {
        warning(
            paste(
                "One of the capture probabilities is zero.",
                "Cannot compute conditional capture probability."
            )
        )
        return(NA)
    }
    gamma_hat <- q_jk / (q_j * q_k)
    return(gamma_hat)
}

#' Get the capture probability estimate given the specified estimator type.
capture_probability_estimate <- function(
    qhat_list,
    estimator = c("PI", "DR", "TMLE")
) {
    estimator <- match.arg(estimator)
    return(switch(
        estimator,
        plugin = get_plugin_estimation(qhat_list),
        double_robust = get_double_robust_estimation(qhat_list),
        tmle = get_tmle_estimation(qhat_list),
        stop(paste(
            "Unsupported estimator specified: ",
            estimator,
            ". Supported options are 'PI', 'DR', and 'TMLE'."
        ))
    ))
}

#' Get the capture probability estimate for the standard Plugin Estimator.
#'
#' @param qhat_list List with elements q_j, q_k, q_jk (numeric vectors)
#' @return The capture probability estimate (psi_hat)
get_plugin_estimation <- function(q_j, q_k, q_jk) {
    gamma_hat <- conditional_capture_prob(q_j, q_k, q_jk)
    capture_prob <- 1 / mean(1 / gamma_hat)
    return(capture_prob)
}


#' Doubly-robust estimator for inverse capture
#' @param qhat_list List with elements q_j, q_k, q_jk (numeric vectors)
#' @param pair Character or vector identifying the pair
#' (not used in calculation)
#' @return List with n_hat, sd_n_hat, psi_inverse_hat, se_psi_inverse_hat
get_double_robust_estimation <- function(q_j, q_k, q_jk, y_j, y_k) {
    n_obs <- length(q_j)

    if (is.null(y_j) || is.null(y_k)) {
        warning(paste0(
            "y_j and y_k not found in qhat_list for pair ",
            pair,
            ". DR estimator reduces to PI estimator."
        ))
        return(get_plugin_estimation(q_j, q_k, q_jk))
    }
    # q_j = observational probability of appearing in list 1
    # q_k = observational probability of appearing in list 2
    # q_jk = observational probability of appearing in both lists
    # y_j = indicator of appearing in list 1
    # y_k = indicator of appearing in list 2
    # gamma = conditional capture probability
    # psi = marginal capture probability
    # psi_hat = estimated marginal capture probability
    # P = theoretical distribution of the entire population (captured and
    #     uncaptured)
    # Q = theoretical distribution of the captured population
    # Q_n = empirical measure under Q: just average the function across the rows

    # DR estimator core calculation
    gamma_hat_inverse <- 1 / conditional_capture_prob(q_j, q_k, q_jk)

    phi_hat <- 1 /
        mean(
            gamma_hat_inverse *
                (y_j / q_j + y_k / q_k - (y_j * y_k) / q_jk)
        )
    return(phi_hat)
}

get_tmle_estimation <- function(qhat_list, pair) {
    return(list(
        n_hat = NA,
        sd_n_hat = NA,
        psi_inverse_hat = NA,
        se_psi_inverse_hat = NA
    ))
}
