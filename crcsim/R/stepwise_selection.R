################################################################################
# File: stepwise_selection.R                                                   #
# Project: crc-sim                                                             #
# Created Date: 2026-05-14                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-06-16                                                    #
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
#' @param model_data a data frame containing the observed capture histories and
#' a frequency column.
#' @param opts a \code{StepwiseOptions} object specifying the options for
#' stepwise regression, including the model family, p-value threshold for
#' variable inclusion, and stepwise direction.
#'
#' @keywords internal
#' @export
stepwise_selection <- function(model_data, opts, verbose = FALSE) {
    if (!inherits(opts, "StepwiseOptions")) {
        stop("Invalid StepwiseOptions object provided.")
    }

    if (!is_frequency_table(model_data, opts$frequency_col_name)) {
        stop(paste(
            "Data must be a frequency table with a numeric frequency column",
            opts$frequency_col_name
        ))
    }
    output <- step_regression(
        model_data,
        opts$frequency_col_name,
        opts$capture_columns,
        p_threshold = opts$threshold,
        direction = opts$direction,
        model_family = opts$model,
        k = opts$interaction_limit,
        verbose = verbose
    )

    return(output)
}


#' Helper function for stepwise regression
#' @param model_data dataframe
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
#' @importFrom purrr quietly
#'
#' @keywords internal
#' @noRd
step_regression <- function(
    model_data,
    y,
    x,
    model_family = c("poisson", "negbin"),
    direction = c("both", "backward", "forward"),
    p_threshold = 0.05,
    k = 2,
    verbose = FALSE
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
    model <- tryCatch(
        fit_loglinear_model(model_data, formula_init, model_family),
        error = function(error) {
            warning(
                sprintf(
                    "Stepwise %s model did not converge: %s",
                    model_family,
                    conditionMessage(error)
                ),
                call. = FALSE
            )
            return(NULL)
        }
    )

    if (is.null(model)) {
        return(stepwise_failure_result(model_family, formula_init))
    }

    # Step function is obnoxious, so we turn off the messages. This requires
    # linter adjustment
    # nolint start: implicit_assignment_linter
    final_model <- tryCatch(
        suppressMessages(
            step(
                model,
                scope = list(upper = formula_max, lower = formula_init),
                direction = direction,
                k = log(nrow(model_data)),
                trace = as.integer(verbose)
            )
        ),
        error = function(error) {
            warning(
                sprintf(
                    "Stepwise %s model selection did not converge: %s",
                    model_family,
                    conditionMessage(error)
                ),
                call. = FALSE
            )
            return(NULL)
        }
    )
    # nolint end: implicit_assignment_linter

    if (is.null(final_model)) {
        return(stepwise_failure_result(model_family, formula_init))
    }

    intercept <- coef(final_model)[1]
    estimate <- exp(intercept)

    # confint forces messages to output which is obnoxious. We suppress them and
    # turn off the linter for the implicit assignment in this block.

    # nolint start: implicit_assignment_linter
    suppressMessages(ci <- exp(confint(final_model)[1, ]))
    # nolint end: implicit_assignment_linter

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

stepwise_failure_result <- function(model_family, formula_init) {
    return(list(
        model = model_family,
        formula = formula_init,
        summary = NULL,
        estimate = NA_real_,
        lower_ci = NA_real_,
        upper_ci = NA_real_,
        AIC = NA_real_
    ))
}
