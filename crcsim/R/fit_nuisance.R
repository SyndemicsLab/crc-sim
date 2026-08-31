################################################################################
# File: fit_nuisance.R                                                         #
# Project: crcsim                                                              #
# Created Date: 2026-08-27                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Estimate nuisance parameters using the specified function.
#'
#' @param func The function to use for nuisance estimation.
#' @param train The training dataset.
#' @param test The testing dataset.
#' @param n_lists The number of lists.
#' @param j Index for the first variable.
#' @param k Index for the second variable.
#' @param margin The margin for estimation.
#'
#' @return Estimated nuisance parameters.
#' @keywords internal
nuisance_estimation <- function(func, train, test, n_lists, j, k, margin) {
    # Temporary implementation only using logit estimation
    return(qhat_logit(train, test, n_lists, j, k, margin))
}


#' Estimate the initial nuisance parameters using GLM fits on the training data
#' and predict on the test data.
#'
#' @param train The training dataset.
#' @param test The testing dataset.
#' @param n_lists The number of lists.
#' @param j Index for the first variable.
#' @param k Index for the second variable.
#' @param margin The margin for estimation.
#'
#' @importFrom stats predict glm binomial
#' @return A list containing the initial estimates for q_1, q_2, and q_12.
#' @keywords internal
qhat_logit <- function(train, test, n_lists, j, k, margin) {
    ## Template Functions for GLM and Q-function estimation
    template_glm <- function(form) {
        return(
            try(glm(
                form,
                family = binomial(link = "logit"),
                data = train[, c(j, k, (n_lists + 1):ncol(train))]
            ))
        )
    }

    template_q_j <- function(d_fit, q_12_offset = 0.0) {
        return(pmin(
            pmax(
                q_12_offset + predict(d_fit, newdata = test, type = "response"),
                margin
            ),
            1
        ))
    }

    ## Core Functionality Start
    c_names <- c(paste0("d", 1:n_lists), paste0("x", 1:(ncol(train) - n_lists)))
    colnames(train) <- c_names
    colnames(test) <- c_names

    fit_j_0 <- template_glm(formula(paste0("d", j, "*(1 - d", k, ") ~.")))
    fit_0_k <- template_glm(formula(paste0("d", k, "*(1 - d", j, ") ~.")))
    fit_j_k <- template_glm(formula(paste0("d", j, "*d", k, " ~.")))

    if (
        inherits(fit_j_0, "try-error") ||
            inherits(fit_0_k, "try-error") ||
            inherits(fit_j_k, "try-error")
    ) {
        warning("One or more GLM fits failed.")
        return(NULL)
    }

    q_12 <- template_q_j(fit_j_k)
    q_1 <- template_q_j(fit_j_0, q_12_offset = q_12)
    q_2 <- template_q_j(fit_0_k, q_12_offset = q_12)

    return(list(q_1 = q_1, q_2 = q_2, q_12 = q_12))
}

#' Fit TMLE nuisance parameters. This updates the initial estimates using the
#' TMLE procedure.
#'
#' @param q_1 Initial estimate for the first nuisance parameter.
#' @param q_2 Initial estimate for the second nuisance parameter.
#' @param q_12 Initial estimate for the joint nuisance parameter.
#' @param y_j Observed outcome for the first variable.
#' @param y_k Observed outcome for the second variable.
#' @param y_jk Observed joint outcome for the first and second variables.
#' @param iterations Maximum number of iterations for the TMLE update.
#' @param margin Margin for estimation.
#' @param n_lists Number of lists.
#'
#' @return Updated nuisance parameter estimates.
#' @keywords internal
tmle_nuisance <- function(
    q_1,
    q_2,
    q_12,
    y_j,
    y_k,
    y_jk,
    iterations,
    margin,
    n_lists
) {
    ############################################################################
    # Helper Functions
    ############################################################################
    logit <- function(x) {
        return(log(x / (1 - x)))
    }

    ############################################################################
    # Reset Nuisance Values
    ############################################################################

    margin_error <- 1 + margin
    count <- 0

    q_10 <- pmin(pmax(q_1 - q_12, margin), 1 - margin)
    q_02 <- pmin(pmax(q_2 - q_12, margin), 1 - margin)
    q_12 <- pmin(pmax(q_12, margin), 1 - margin)

    y_j0 <- y_j * (1 - y_k)
    y_0k <- (1 - y_j) * y_k

    while (abs(margin_error) > margin && count < iterations) {
        q12_fit <- fit_glm_tmle(
            y_jk,
            logit(q_12),
            (q_10 + q_12) /
                q_12 +
                (q_02 + q_12) / q_12 -
                (q_10 + q_12) * (q_02 + q_12) / q_12^2
        )
        if (!is.null(q12_fit$value)) {
            q_12 <- pmax(pmin(q12_fit$value, 1), margin)
        }

        q10_fit <- fit_glm_tmle(
            y_j0,
            logit(q_10),
            (q_02 + q_12) / q_12
        )
        if (!is.null(q10_fit$value)) {
            q_10 <- pmax(pmin(q10_fit$value, 1 - q_12), margin)
        }

        if (n_lists > 2) {
            q02_fit <- fit_glm_tmle(
                y_0k,
                logit(q_02),
                (q_10 + q_12) / q_12
            )
            if (!is.null(q02_fit$value)) {
                q_02 <- pmax(
                    pmin(q02_fit$value, 1 - q_10 - q_12),
                    margin
                )
            }
        } else {
            q02_fit <- list(error = 0)
            q_02 <- pmax(0, 1 - q_10 - q_12)
        }

        margin_error <- max(q12_fit$error, q10_fit$error, q02_fit$error)
        count <- count + 1
    }

    new_nuisances <- list(
        q_1 = q_10 + q_12,
        q_2 = q_02 + q_12,
        q_12 = q_12
    )
    return(new_nuisances)
}

#' Fit a GLM for TMLE updates.
#'
#' @param response The response variable.
#' @param offset The offset for the GLM.
#' @param ratio The ratio used as the predictor in the GLM.
#'
#' @return A list containing the fitted values and the maximum absolute
#' coefficient value as the error.
#'
#' @importFrom stats complete.cases glm.fit coef binomial
#' @keywords internal
fit_glm_tmle <- function(response, offset, ratio) {
    complete <- complete.cases(response, offset, ratio)

    fit <- try(
        glm.fit(
            x = matrix(ratio[complete], ncol = 1),
            y = response[complete],
            offset = offset[complete],
            family = binomial(link = "logit")
        ),
        silent = TRUE
    )

    if (inherits(fit, "try-error")) {
        return(list(value = NULL, error = Inf))
    }

    return(list(
        value = fit$fitted.values,
        error = max(abs(coef(fit)), na.rm = TRUE)
    ))
}
