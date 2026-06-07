################################################################################
# File: capture_probability.R                                                  #
# Project: crc-sim                                                             #
# Created Date: 2026-05-18                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Estimate Capture Probability
#' @description Estimate the capture probability using the specified estimator
#' (e.g., Plugin, Doubly-Robust, TMLE). This function serves as a wrapper that
#' calls the appropriate estimation function based on the user's choice of
#' estimator. The function takes in the estimated capture probabilities (q_j,
#' q_k, q_jk) and optionally the outcome indicators (y_j, y_k) for the
#' doubly-robust and TMLE estimators. The function checks the validity of the
#' input parameters and returns the estimated capture probability based on the
#' chosen estimator. If an unsupported estimator is specified, an error is
#' raised.
#'
#' @param q_j Numeric vector of capture probabilities for list j
#' @param q_k Numeric vector of capture probabilities for list k
#' @param q_jk Numeric vector of capture probabilities for both lists j and k
#' @param y_j Optional numeric vector of outcome indicators for list j
#' (required for DR and TMLE estimators)
#' @param y_k Optional numeric vector of outcome indicators for list k
#' (required for DR and TMLE estimators)
#' @param estimator Character string specifying the estimator to use ("PI",
#' "DR", or "TMLE")
#' @return The estimated capture probability (psi_hat) based on the chosen
#' estimator
#'
#' @export
estimate_capture_probability <- function(
    q_j,
    q_k,
    q_jk,
    y_j = NULL,
    y_k = NULL,
    estimator = c("PI", "DR", "TMLE")
) {
    estimator <- match.arg(estimator)
    return(switch(
        estimator,
        PI = get_plugin_estimation(q_j, q_k, q_jk),
        DR = get_double_robust_estimation(q_j, q_k, q_jk, y_j, y_k),
        TMLE = get_tmle_estimation(q_j, q_k, q_jk, y_j, y_k),
        stop(paste(
            "Unsupported estimator specified: ",
            estimator,
            ". Supported options are 'PI', 'DR', and 'TMLE'."
        ))
    ))
}

#' Conditional Capture Probability Calculation
#' @description Calculate the conditional capture probability (gamma_hat) given
#' q_j, q_k, and q_jk. If we are able to assume that q_j is the capture
#' probability when Y_j == 1 conditional on covariate set X, q_k is Y_k == 1
#' given X, and q_jk is Y_j == 1 and Y_k == 1 given X, then the conditional
#' capture probability is the probability that Y != 0 given X. This function is
#' used in the estimation of capture probabilities using various estimators
#' (e.g., Plugin, Doubly-Robust). If either q_j or q_k is zero, a warning is
#' issued and NA is returned since the conditional capture probability cannot
#' be computed. The function assumes that q_j, q_k, and q_jk are of the same
#' length and correspond to the same set of observations.
#'
#' @param q_j Numeric vector of capture probabilities for list j
#' @param q_k Numeric vector of capture probabilities for list k
#' @param q_jk Numeric vector of capture probabilities for both lists j and k
#' @return Numeric vector of conditional capture probabilities (gamma_hat)
#'
#' @keywords internal
#' @export
conditional_capture_prob <- function(q_j, q_k, q_jk) {
    if (
        is_probability_vector(q_j) &&
            is_probability_vector(q_k) &&
            is_probability_vector(q_jk)
    ) {
        if (length(q_j) != length(q_k) || length(q_j) != length(q_jk)) {
            warning("q_j, q_k, and q_jk must be of the same length.")
            return(NA)
        }
    } else if (
        !is_probability_scalar(q_j) ||
            !is_probability_scalar(q_k) ||
            !is_probability_scalar(q_jk)
    ) {
        warning("q_j, q_k, and q_jk must all be between 0 and 1.")
        return(NA)
    }

    gamma_hat <- q_jk / (q_j * q_k)
    return(gamma_hat)
}

#' Plugin estimator for capture probability
#' @description Get the capture probability estimate for the standard plugin
#' estimator. This estimator calculates the conditional capture probability
#' (gamma_hat) for each observation and then computes the harmoinc mean of the
#' inverse of gamma_hat to get the final capture probability estimate. The
#' plugin estimator relies solely on the estimated capture probabilities
#' (q_j, q_k, q_jk) without incorporating any outcome information (y_j, y_k).
#'
#' @param q_j Numeric vector of capture probabilities for list j
#' @param q_k Numeric vector of capture probabilities for list k
#' @param q_jk Numeric vector of capture probabilities for both lists j and k
#' @return The capture probability estimate (psi_hat)
#'
#' @keywords internal
#' @export
get_plugin_estimation <- function(q_j, q_k, q_jk) {
    gamma_hat <- conditional_capture_prob(q_j, q_k, q_jk)
    capture_prob <- 1 / mean(1 / gamma_hat)
    return(capture_prob)
}


#' Doubly-robust estimator for capture probability
#' @description Get the capture probability estimate for the doubly-robust
#' estimator. This estimator is built off the plugin estimator idea but
#' incorporates outcome information (y_j, y_k) to debias the estimate. The
#' debiasing is done by calculating an efficient influence function (phi_hat).
#' For information about the derivation of the doubly-robust estimator, see the
#' following paper: https://arxiv.org/abs/2104.14091. If either y_j or y_k is
#' not provided, a warning is issued and the function falls back to the plugin
#' estimator since the doubly-robust estimator cannot be computed without
#' outcome information.
#'
#' @param qhat_list List with elements q_j, q_k, q_jk (numeric vectors)
#' @param pair Character or vector identifying the pair
#' (not used in calculation)
#' @return List with n_hat, sd_n_hat, psi_inverse_hat, se_psi_inverse_hat
#'
#' @references Das, M., Kennedy, E. H., & Jewell, N.P. (2021). Doubly robust
#' capture-recapture methods for estimating population size. _arXiv preprint_
#' arXiv:2104.14091.
#'
#' @keywords internal
#' @export
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
    gamma_hat_inverse <- 1 / conditional_capture_prob(q_j, q_k, q_jk)

    capture_prob <- 1 /
        mean(
            gamma_hat_inverse *
                (y_j / q_j + y_k / q_k - (y_j * y_k) / q_jk)
        )
    return(capture_prob)
}


#' TMLE estimator for capture probability
#' @description Get the capture probability estimate for the TMLE estimator.
#' This estimator builds off the doubly-robust estimator by iteratively
#' updating the estimates of q_j, q_k, and q_jk using the clever covariate
#' until convergence.
#'
#' @param q_init_j Numeric vector of initial capture probabilities for list j
#' @param q_init_k Numeric vector of initial capture probabilities for list k
#' @param q_init_jk Numeric vector of initial capture probabilities for both
#' lists j and k
#' @param y_j Numeric vector of outcome indicators for list j
#' @param y_k Numeric vector of outcome indicators for list k
#' @param iter Integer specifying the maximum number of iterations for fitting
#' @param error Numeric value specifying the convergence threshold for the
#' clever covariate fitting process
#' @return The estimated capture probability (psi_hat) based on the TMLE
#' estimator
#'
#' @references van der Laan, M. J. and Rubin, D. (2006). Targeted maximum
#' likelihood learning. _The International Journal of Biostatistics_, *2*(1)
#' @references Das, M., Kennedy, E. H., & Jewell, N.P. (2021). Doubly robust
#' capture-recapture methods for estimating population size. _arXiv preprint_
#' arXiv:2104.14091.
#'
#' @keywords internal
#' @export
get_tmle_estimation <- function(
    q_init_j,
    q_init_k,
    q_init_jk,
    y_j,
    y_k,
    iter = 250,
    error = 0.005
) {
    # Placeholder for TMLE estimation logic
    clever_covariate <- fit_clever_covariate(
        q_init_j,
        q_init_k,
        q_init_jk,
        y_j,
        y_k,
        iter,
        error
    )
    q_conv_j <- clever_covariate[["q_conv_j"]]
    q_conv_k <- clever_covariate[["q_conv_k"]]
    q_conv_jk <- clever_covariate[["q_conv_jk"]]
    capture_prob <- 1 / mean((q_conv_j * q_conv_k) / q_conv_jk)
    return(capture_prob)
}


#' Expit function
#' @description The expit function, also known as the logistic function, is used
#' to transform log-odds (logit) back to probabilities. It is defined as:
#' exp(x) / (1 + exp(x)).
#'
#' @param x Numeric vector of log-odds values
#' @return Numeric vector of probabilities corresponding to the input log-odds
#'
#' @keywords internal
#' @noRd
expit <- function(x) {
    return(exp(x) / (1 + exp(x)))
}

#' Logit function
#' @description The logit function is the inverse of the expit function and is
#'  used to transform probabilities to log-odds. It is defined as:
#' log(p / (1 - p)).
#'
#' @param p Numeric vector of probabilities (values between 0 and 1)
#' @return Numeric vector of log-odds corresponding to the input probabilities
#'
#' @keywords internal
#' @noRd
logit <- function(p) {
    if (any(p >= 1)) {
        stop("Input cannot be 1 since logit is undefined for p = 1.")
    }
    return(log(p / (1 - p)))
}

prob_bound <- function(p) {
    return(pmin(pmax(p, 1e-6), 1 - 1e-6))
}

#' Fit clever covariate for TMLE estimation
#' @description Fit the clever covariate for TMLE estimation of capture
#' probability. This function iteratively updates the estimates of q_j, q_k,
#' and q_jk using the clever covariates until convergence or until the maximum
#' number of iterations is reached.
#'
#' @references van der Laan, M. J. and Rubin, D. (2006). Targeted maximum
#' likelihood learning. _The International Journal of Biostatistics_, *2*(1)
#' @references Das, M., Kennedy, E. H., & Jewell, N.P. (2021). Doubly robust
#' capture-recapture methods for estimating population size. _arXiv preprint_
#' arXiv:2104.14091.
#'
#' @param q_init_j Numeric vector of initial capture probabilities for list j
#' @param q_init_k Numeric vector of initial capture probabilities for list k
#' @param q_init_jk Numeric vector of initial capture probabilities for both
#' lists j and k
#' @param y_j Numeric vector of outcome indicators for list j
#' @param y_k Numeric vector of outcome indicators for list k
#' @param iter Integer specifying the maximum number of iterations for fitting
#' the clever covariate
#' @param error Numeric value specifying the convergence threshold for the
#' clever covariate fitting process
#' @return A list containing the converged estimates of q_j, q_k, and q_jk
#' after fitting the clever covariate
#'
#' @keywords internal
fit_clever_covariate <- function(
    q_init_j,
    q_init_k,
    q_init_jk,
    y_j,
    y_k,
    iter,
    error
) {
    q_j <- q_init_j
    q_k <- q_init_k
    q_jk <- q_init_jk
    counter <- 0
    converged <- FALSE

    y_j0 <- y_j * (1 - y_k)
    y_k0 <- y_k * (1 - y_j)

    while (!converged && counter < iter) {
        h_j <- (q_k / q_jk)
        h_k <- (q_j / q_jk)
        h_jk <- ((q_j * q_k) / q_jk^2) - h_k - h_j

        # Now we regress 3 models using the clever covariates as predictors

        # Model 1 models the next step for q_jk
        m_1 <- glm(
            y_j ~ offset(logit(q_j)) + h_jk - 1,
            family = binomial(link = "logit"),
            data = data.frame(
                y_j = y_j,
                q_j = q_j,
                h_jk = h_jk
            ),
            na.action = na.omit
        )
        q_next_jk <- prob_bound(expit(logit(q_jk) + coef(m_1) * h_jk))
        if (any(q_next_jk <= 0 | q_next_jk >= 1)) {
            stop(
                paste(
                    "Estimated q_jk values must be between 0 and 1.",
                    "Check model fit. Estimated q_jk values are:",
                    q_next_jk
                )
            )
        }

        # Model 2 models the next step for q_j
        m_2 <- glm(
            y_j0 ~ offset(logit(q_j - q_next_jk)) + h_j - 1,
            family = binomial(link = "logit"),
            data = data.frame(
                y_j0 = y_j0,
                q_j = q_j,
                h_j = h_j
            ),
            na.action = na.omit
        )
        q_next_j <- prob_bound(pmin(
            q_next_jk + expit(logit(q_j - q_next_jk) + coef(m_2) * h_j),
            1 - q_next_jk
        ))
        if (any(q_next_j <= 0 | q_next_j >= 1)) {
            stop(
                paste(
                    "Estimated q_j values must be between 0 and 1.",
                    "Check model fit. Estimated q_j values are:",
                    q_next_j
                )
            )
        }

        m_3 <- glm(
            y_k0 ~ offset(logit(q_k - q_next_jk)) + h_k - 1,
            family = binomial(link = "logit"),
            data = data.frame(
                y_k0 = y_k0,
                q_k = q_k,
                h_k = h_k
            ),
            na.action = na.omit
        )
        q_next_k <- prob_bound(pmin(
            q_next_jk + expit(logit(q_k - q_next_jk) + coef(m_3) * h_k),
            1 + q_next_jk - q_next_j
        ))
        if (any(q_next_k <= 0 | q_next_k >= 1)) {
            stop(
                paste(
                    "Estimated q_k values must be between 0 and 1.",
                    "Check model fit. Estimated q_k values are:",
                    q_next_k
                )
            )
        }

        q_j <- q_next_j
        q_k <- q_next_k
        q_jk <- q_next_jk

        if (max(coef(m_1), coef(m_2), coef(m_3)) <= error) {
            converged <- TRUE
        }
        counter <- counter + 1
    }
    # Placeholder for clever covariate calculation logic
    return(list(
        q_conv_j = q_j,
        q_conv_k = q_k,
        q_conv_jk = q_jk
    ))
}
