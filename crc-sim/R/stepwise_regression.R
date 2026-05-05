################################################################################
# File: stepwise_regression.R                                                  #
# Project: crc-sim                                                             #
# Created Date: 2026-05-05                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-05                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Helper function for stepwise regression
#' @param data dataframe
#' @param y string: LHS of formula object
#' @param x string: RHS of formula object
#' @param method string: either 'poisson' or 'negbin'
#' @param direction string: stepwise direction
#' @param p_threshold numeric: threshold for stepwise selection
#' @param k integer: limit for k-way interaction terms
#' @param verbose logical: whether to print intermediate models
#'
#' @keywords internal
#' @importFrom MASS glm.nb
#' @importFrom utils capture.output
#' @importFrom stats AIC coef confint formula glm poisson step

step_regression <- function(
    data,
    y,
    x,
    method = "poisson",
    direction = "both",
    p_threshold = 0.05,
    k = 2,
    verbose = TRUE
) {
    formula_init <- as.formula(paste(y, "~", paste(x, collapse = " + ")))
    formula_max <- as.formula(paste(
        y,
        "~ (",
        paste(x, collapse = " + "),
        ")^",
        k
    ))

    if (verbose) {
        if (method == "poisson") {
            init_mod <- glm(formula_init, family = poisson, data = data)
        } else {
            init_mod <- MASS::glm.nb(formula_init, data = data)
        }

        final_mod <- step(
            init_mod,
            scope = list(upper = formula_max, lower = formula_init),
            direction = direction,
            k = log(nrow(data))
        )
    } else {
        capture.output({
            if (method == "poisson") {
                init_mod <- glm(formula_init, family = poisson, data = data)
            } else {
                init_mod <- MASS::glm.nb(formula_init, data = data)
            }

            final_mod <- suppressWarnings(step(
                init_mod,
                scope = list(upper = formula_max, lower = formula_init),
                direction = direction,
                k = log(nrow(data))
            ))
        })
    }

    intercept <- coef(final_mod)[1]
    estimate <- exp(intercept)
    ci <- exp(confint(final_mod)[1, ])

    results <- list(
        model = method,
        formula = formula(final_mod),
        summary = summary(final_mod),
        estimate = unname(round(estimate, 2)),
        lower_ci = unname(round(ci[1], 2)),
        upper_ci = unname(round(ci[2], 2)),
        AIC = AIC(final_mod)
    )

    return(results)
}

fit_model <- function(
    method,
    formula_init,
    formula_max,
    data,
    direction
) {
    if (method == "poisson") {
        init_mod <- glm(formula_init, family = poisson, data = data)
    } else {
        init_mod <- MASS::glm.nb(formula_init, data = data)
    }

    return(step(
        init_mod,
        scope = list(upper = formula_max, lower = formula_init),
        direction = direction,
        k = log(nrow(data))
    ))
}
