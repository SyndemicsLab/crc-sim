################################################################################
# File: stepwise_selection.R                                                   #
# Project: crc-sim                                                             #
# Created Date: 2026-05-14                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-14                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Run CRC with Stepwise Formula Selection.
#'
#' @keywords internal
#' @export
run_crc_stepwise <- function(
    data,
    freq_column,
    binary_variables,
    method,
    opts_stepwise
) {
    return(step_regression(
        data,
        freq_column,
        binary_variables,
        p_threshold = opts_stepwise$threshold,
        direction = opts_stepwise$direction,
        method = method,
        verbose = opts_stepwise$verbose,
        k = 2
    ))
}


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
#' @importFrom MASS glm.nb
#' @importFrom utils capture.output
#' @importFrom stats AIC coef confint formula glm poisson step
#'
#' @keywords internal
#' @noRd
step_regression <- function(
    data,
    y,
    x,
    model_family = c("poisson", "negbin"),
    direction = "both",
    p_threshold = 0.05,
    k = 2,
    verbose = TRUE
) {
    model_family <- match.arg(model_family)
    formula_init <- as.formula(paste(y, "~", paste(x, collapse = " + ")))
    formula_max <- as.formula(paste(
        y,
        "~ (",
        paste(x, collapse = " + "),
        ")^",
        k
    ))

    model <- fit_loglinear_model(data, formula_init, model_family)

    step(
        model,
        scope = list(upper = formula_max, lower = formula_init),
        direction = direction,
        k = log(nrow(data))
    )

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
