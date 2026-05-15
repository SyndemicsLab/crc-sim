################################################################################
# File: crc_stepwise.R                                                         #
# Project: crc-sim                                                             #
# Created Date: 2026-05-14                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Run CRC with Stepwise Formula Selection.
#' @description This function runs CRC using stepwise regression for formula
#' selection. It takes in the data and options for stepwise regression, fits the
#' specified model using stepwise regression, and returns a list with the final
#' formula, estimate, confidence interval, and AIC value for the selected model.
#'
#' @param data a data frame containing the observed capture histories and a
#' frequency column.
#' @param opts a \code{StepwiseOptions} object specifying the options for
#' stepwise regression, including the model family, p-value threshold for
#' variable inclusion, and stepwise direction.
#'
#' @keywords internal
#' @export
crc_stepwise <- function(data, opts) {
    if (!inherits(opts, "StepwiseOptions")) {
        stop("Invalid StepwiseOptions object provided.")
    }

    # TODO: validate the data is a frequency table
    output <- step_regression(
        data,
        opts$freq_column,
        opts$capture_indicators,
        p_threshold = opts$threshold,
        direction = opts$direction,
        method = opts$model,
        k = opts$interaction_limit
    )

    return(output)
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
    direction = c("both", "backward", "forward"),
    p_threshold = 0.05,
    k = 2
) {
    model_family <- match.arg(model_family)
    direction <- match.arg(direction)

    formula_init <- as.formula(paste(y, "~", paste(x, collapse = " + ")))
    formula_max <- as.formula(paste(
        y,
        "~ (",
        paste(x, collapse = " + "),
        ")^",
        k
    ))

    model <- fit_loglinear_model(data, formula_init, model_family)

    final_model <- step(
        model,
        scope = list(upper = formula_max, lower = formula_init),
        direction = direction,
        k = log(nrow(data))
    )

    intercept <- coef(final_model)[1]
    estimate <- exp(intercept)
    ci <- exp(confint(final_model)[1, ])

    results <- list(
        model = model_family,
        formula = formula(final_model),
        summary = summary(final_model),
        estimate = unname(round(estimate, 2)),
        lower_ci = unname(round(ci[1], 2)),
        upper_ci = unname(round(ci[2], 2)),
        AIC = AIC(final_model)
    )

    return(results)
}
