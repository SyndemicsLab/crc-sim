################################################################################
# File: aic_selection.R                                                        #
# Project: crc-sim                                                             #
# Created Date: 2026-05-14                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-21                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Run CRC with AIC Formula Selection. This function takes in the data and a
#' list of potential formulas to fit, fits the specified model for each
#' formula, and returns a dataframe with the formula, estimate, AIC value,
#' confidence interval, and any error messages for each formula.
#'
#' @param data a data frame containing the observed capture histories and a
#' frequency column.
#' @param formulas a list of formula objects specifying the log-linear models to
#' fit and compare using AIC
#' @param model_family a string specifying the model family to use, either
#' "poisson" or "negbin". Default is "poisson".
#' @return a data frame with the formula, estimate, AIC value, confidence
#' interval, and any error messages for each formula
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows arrange
#'
#' @keywords internal
#' @export
aic_selection <- function(data, opts) {
    if (!inherits(opts, "AICOptions")) {
        stop("Invalid AICOptions object provided.")
    }

    if (!is_frequency_table(data, opts$frequency_col_name)) {
        stop(paste(
            "Data must be a frequency table with a numeric frequency column",
            opts$frequency_col_name
        ))
    }

    output <- map(
        opts$formulas,
        evaluate_formula_with_aic,
        data = data,
        model_family = opts$model
    ) |>
        bind_rows() |>
        arrange(.data[["AIC"]])

    return(output)
}

#' Evaluate a single possible formula based on AIC. This function is an
#' internal function called specifically by \code{run_aic_selection} to
#' evaluate each formula in the provided list of formulas. It fits the
#' specified model for the given formula and returns a dataframe with the
#' formula, estimate, AIC value, confidence interval, and any error messages
#' for the formula.
#'
#' @param formula_object a formula object specifying the log-linear model to fit
#' @param data a data frame containing the observed capture histories and a
#' frequency column.
#' @param model_family a string specifying the model family to use, either
#' "poisson" or "negbin". Default is "poisson".
#' @return a data frame with the formula, estimate, AIC value, confidence
#' interval, and any error messages for the specified formula
#'
#' @keywords internal
#' @noRd
evaluate_formula_with_aic <- function(
    formula_object,
    data,
    model_family = c("poisson", "negbin")
) {
    model_family <- match.arg(model_family)
    return(tryCatch(
        {
            model <- fit_loglinear_model(data, formula_object, model_family)
            return(build_aic_success_row(model, formula_object))
        },
        error = function(error) {
            return(build_aic_error_row(formula_object, error))
        }
    ))
}

#' Build AIC Result Row on Success. An internal function used by \code
#' {evaluate_formula_with_aic} to build a result row for the AIC evaluation
#' when the model fitting is successful.
#'
#' @param model a fitted model object resulting from the log-linear model fit
#' @param formula_object a formula object specifying the log-linear model that
#' was fit
#' @return a data frame with the formula, estimate, AIC value, confidence
#' interval, and NA for the error message
#'
#' @importFrom stats coef AIC confint
#' @keywords internal
#' @noRd
build_aic_success_row <- function(model, formula_object) {
    intercept <- exp(coef(model)["(Intercept)"])
    aic_value <- AIC(model)
    ci <- suppressMessages(confint(
        model,
        "(Intercept)",
        level = 0.95
    ))

    return(data.frame(
        formula = paste(deparse(formula_object), collapse = " "),
        estimate = round(intercept, 2),
        AIC = round(aic_value, 2),
        lower_ci = unname(round(exp(ci[1]), 2)),
        upper_ci = unname(round(exp(ci[2]), 2)),
        error = NA,
        row.names = NULL
    ))
}

#' Build AIC Result Row on Error. An internal function used by
#' \code{evaluate_formula_with_aic} to build a result row for the AIC
#' evaluation when there is an error in model fitting.
#'
#' @param formula_object a formula object specifying the log-linear model that
#' was attempted to be fit
#' @param error an error object containing the error message from the failed
#' model fitting attempt
#' @return a data frame with the formula, NA for the estimate, AIC value, and
#' confidence interval, and the error message
#'
#' @keywords internal
#' @noRd
build_aic_error_row <- function(formula_object, error) {
    return(data.frame(
        formula = paste(deparse(formula_object), collapse = " "),
        estimate = NA,
        AIC = NA,
        lower_ci = NA,
        upper_ci = NA,
        error = toString(error$message),
        row.names = NULL
    ))
}
