################################################################################
# File: crc.R                                                                  #
# Project: crc-sim                                                             #
# Created Date: 2026-02-23                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-02-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Spatial Capture Re-Capture
#'
#' A method for estimation of 'unknowns' through knowledge about knowns as
#' described by Barocas, Joshua A et al. "Estimated Prevalence of Opioid Use
#' Disorder in Massachusetts, 2011-2015: A Capture-Recapture Analysis."
#' doi:10.2105/AJPH.2018.304673
#'
#' This implementation supports poisson and negative binomial regression models
#' with either AIC-based formula selection or stepwise selection.
#'
#' @param data Dataframe: A dataframe containing a frequency column and binary
#' columns indicating involvement in the given database
#' @param freq.column Column: A column containing the frequency of observed
#' combinations
#' @param binary.variables List of Columns: List containing columns of binary
#' variables indicating involvement in the given database
#' @param method String: Selection for the spatial capture-recapture method -
#' either 'poisson' or 'negbin'
#' @param formula.selection String: Selection for formula decision - either
#' 'aic' or 'stepwise'
#' @param formula Formula: Optional custom formula object for regression
#' @param opts.stepwise List: List of \code{direction}: 'forward', 'backward',
#' or 'both', \code{threshold}: p-value threshold for stepwise selection, and
#' \code{verbose} if you would like every stepped-through model to be printed
#'
#' @import data.table
#' @importFrom MASS glm.nb
#' @importFrom stats AIC coef confint formula glm poisson
#' @export

crc <- function(
    data,
    freq.column,
    binary.variables,
    method = "poisson",
    formula.selection = "stepwise",
    formula = NULL,
    opts.stepwise = list(
        direction = "both",
        threshold = 0.05,
        verbose = TRUE
    )
) {
    if (!(method %in% c("poisson", "negbin"))) {
        stop("Method must be either 'poisson' or 'negbin'")
    }

    if (!(formula.selection %in% c("aic", "stepwise"))) {
        stop("Formula selection must be either 'aic' or 'stepwise'")
    }

    dt <- setDT(data)

    if (is.null(formula)) {
        if (formula.selection == "aic") {
            form <- formula_list(freq.column, binary.variables)
        } else if (formula.selection == "stepwise") {
            form <- NULL
        }
    } else if (!is.formula(formula)) {
        stop("Expected Formula Object when Specifying Formula")
    } else {
        form <- formula
    }

    if (formula.selection == "aic") {
        results <- list()

        for (i in seq_along(form)) {
            result <- tryCatch(
                {
                    if (method == "poisson") {
                        model <- stats::glm(
                            form[[i]],
                            data = dt,
                            family = "poisson"
                        )
                    } else {
                        model <- MASS::glm.nb(formula = form[[i]], data = dt)
                    }

                    intercept <- exp(stats::coef(model)["(Intercept)"])
                    aic <- stats::AIC(model)

                    ci <- suppressMessages(stats::confint(
                        model,
                        "(Intercept)",
                        level = 0.95
                    ))
                    lower_ci <- exp(ci[1])
                    upper_ci <- exp(ci[2])

                    data.frame(
                        formula = paste(deparse(form[[i]]), collapse = " "),
                        estimate = round(intercept, 2),
                        AIC = round(aic, 2),
                        lower_ci = unname(round(lower_ci, 2)),
                        upper_ci = unname(round(upper_ci, 2)),
                        error = NA,
                        row.names = NULL
                    )
                },
                error = function(e) {
                    data.frame(
                        formula = paste(deparse(form[[i]]), collapse = " "),
                        estimate = NA,
                        AIC = NA,
                        lower_ci = NA,
                        upper_ci = NA,
                        error = toString(e$message),
                        row.names = NULL
                    )
                }
            )

            results[[i]] <- result
        }

        model <- do.call(rbind, results)
        model <- model[order(model$AIC), ]
    } else if (formula.selection == "stepwise") {
        model <- step_regression(
            data,
            freq.column,
            binary.variables,
            p.threshold = opts.stepwise$threshold,
            direction = opts.stepwise$direction,
            method = method,
            verbose = opts.stepwise$verbose,
            k = 2
        )
    }

    return(model)
}
