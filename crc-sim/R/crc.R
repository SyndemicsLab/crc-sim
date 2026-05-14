################################################################################
# File: crc.R                                                                  #
# Project: crc-sim                                                             #
# Created Date: 2026-02-23                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-14                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Validate CRC Inputs
#'
#' @keywords internal
#' @noRd

validate_crc_args <- function(method, formula_selection, formula) {
    if (!(method %in% c("poisson", "negbin"))) {
        stop("Method must be either 'poisson' or 'negbin'")
    }

    if (!(formula_selection %in% c("aic", "stepwise", "tmle"))) {
        stop("Formula selection must be one of: 'aic', 'stepwise', 'tmle'")
    }

    if (!is.null(formula) && !is_formula(formula)) {
        stop("Expected Formula Object when Specifying Formula")
    }

    return(NULL)
}

#' Resolve Formula Input for CRC
#'
#' @keywords internal
#' @noRd

resolve_formula_input <- function(
    freq_column,
    binary_variables,
    formula_selection,
    formula
) {
    if (!is.null(formula)) {
        return(formula)
    }

    if (formula_selection == "aic") {
        return(formula_list(freq_column, binary_variables))
    }

    return(NULL)
}


#' Run CRC with Targeted Maximum Likelihood Estimation (TMLE)
#'
#' @keywords internal
#' @noRd

run_crc_tmle_wrapper <- function(
    data,
    freq_column,
    binary_variables,
    opts_tmle
) {
    return(run_crc_tmle(
        data = data,
        freq_column = freq_column,
        binary_variables = binary_variables,
        opts_tmle = opts_tmle
    ))
}

#' Spatial Capture Re-Capture
#'
#' A method for estimation of 'unknowns' through knowledge about knowns as
#' described by Barocas, Joshua A et al. "Estimated Prevalence of Opioid Use
#' Disorder in Massachusetts, 2011-2015: A Capture-Recapture Analysis."
#' doi:10.2105/AJPH.2018.304673
#'
#' This implementation supports poisson and negative binomial regression models
#' with either AIC-based, stepwise, or TMLE formula selection.
#'
#' @param data Dataframe: A dataframe containing a frequency column and binary
#' columns indicating involvement in the given database
#' @param freq_column Column: A column containing the frequency of observed
#' combinations
#' @param binary_variables List of Columns: List containing columns of binary
#' variables indicating involvement in the given database
#' @param method String: Selection for the spatial capture-recapture method -
#' either 'poisson' or 'negbin'
#' @param formula_selection String: Selection for formula decision - either
#' 'aic', 'stepwise', or 'tmle'
#' @param formula Formula: Optional custom formula object for regression
#' @param opts_stepwise List: List of \code{direction}: 'forward', 'backward',
#' or 'both', \code{threshold}: p-value threshold for stepwise selection, and
#' \code{verbose} if you would like every stepped-through model to be printed
#' @param opts_tmle List: TMLE options where \code{list_pair} is required and
#' should specify two capture columns (e.g., \code{c("capture_1","capture_2")}).
#' Optional fields include \code{funcname}, \code{nfolds}, \code{margin},
#' \code{estimator}, \code{expansion_mode}, \code{sample_size}, and
#' \code{warn_expanded_rows}
#'
#' @importFrom MASS glm.nb
#' @importFrom drpop popsize
#' @importFrom stats AIC coef confint formula glm poisson
#' @export

crc <- function(
    data,
    freq_column,
    binary_variables,
    method = "poisson",
    formula_selection = "stepwise",
    formula = NULL,
    opts_stepwise = list(
        direction = "both",
        threshold = 0.05,
        verbose = TRUE
    ),
    opts_tmle = list(
        list_pair = NULL,
        funcname = "logit",
        nfolds = 2,
        margin = 0.005,
        estimator = "TMLE",
        expansion_mode = "exact",
        sample_size = 100000,
        warn_expanded_rows = 1000000
    )
) {
    validate_crc_args(method, formula_selection, formula)
    data <- tibble::as_tibble(data)
    formulas <- resolve_formula_input(
        freq_column,
        binary_variables,
        formula_selection,
        formula
    )

    if (formula_selection == "aic") {
        return(run_crc_aic(data, formulas, method))
    }

    if (formula_selection == "tmle") {
        return(run_crc_tmle_wrapper(
            data,
            freq_column,
            binary_variables,
            opts_tmle
        ))
    }

    return(run_crc_stepwise(
        data,
        freq_column,
        binary_variables,
        method,
        opts_stepwise
    ))
}
