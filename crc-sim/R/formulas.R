################################################################################
# File: formulas.R                                                             #
# Project: crc-sim                                                             #
# Created Date: 2026-05-05                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-21                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Is Formula
#' @description Generic test of an object being interpretable as a formula
#'
#' @param x object to be tested
#' @return TRUE if object inherits from class "formula", FALSE otherwise
#'
#' @export
is_formula <- function(x) {
    return(inherits(x, "formula"))
}

#' Generates All Possible Combination of Interaction Terms
#' Function that returns a list of every possible combination of
#' interaction terms
#'
#' @param y Column: Column on the estimation side of the equation.
#' @param x List: Columns to return every combination of
#' @importFrom utils combn
#' @importFrom stats as.formula
#'
#' @export
# nolint start: cyclocomplexity_linter
formula_list <- function(y, x) {
    n <- length(x)
    all_formulas <- list()

    for (i in 1:n) {
        all_formulas <- c(all_formulas, paste0(y, "~", x[i]))
    }

    for (i in 2:n) {
        combinations <- combn(x, i)
        for (j in seq_len(ncol(combinations))) {
            combination <- combinations[, j]

            all_formulas <- c(
                all_formulas,
                paste0(y, "~", paste(combination, collapse = "+"))
            )
            interaction_formula <- paste(combination, collapse = "*")
            all_formulas <- c(all_formulas, paste0(y, "~", interaction_formula))

            for (k in 1:(i - 1)) {
                combinations_additive <- combn(combination, k)
                for (l in seq_len(ncol(combinations_additive))) {
                    combination_additive <- combinations_additive[, l]
                    all_formulas <- c(
                        all_formulas,
                        paste0(
                            y,
                            "~",
                            paste(combination_additive, collapse = "+"),
                            "+",
                            interaction_formula
                        )
                    )
                }
            }
        }
    }

    return(lapply(unique(all_formulas), as.formula))
}
# nolint end: cyclocomplexity_linter
