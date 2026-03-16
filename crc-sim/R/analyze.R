################################################################################
# File: analyze.R                                                              #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-02-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Function to run the analysis
#'
#' @param data_table data.table: data.table created from \code{create_data}
#' @param suppress int: supression cap, defaults to 10
#' @param method string: selection for spatial capture recapture model
#'  - either "poisson" or "negbin"
#' @param formula_selection string: selection for formula decision
#'  - either "aic", "corr", or "stepwise"
#' @param opts_stepwise list: list containing 'direction' of 'forward'
#'  'backward' or 'both', and 'threshold': p value threshold for stepwise
#'  selection
#'
#' @import data.table
#' @keywords internal
analyze <- function(
    data_table,
    suppress,
    method = "poisson",
    formula_selection = "stepwise",
    opts_stepwise = list(
        direction = "forward",
        threshold = 0.5,
        verbose = FALSE
    )
) {
    captures <- names(data_table)[
        !names(data_table) %in% c("N_ID", "strata", "tmp")
    ]

    ground_truth <- extract_ground_truth(
        data_table,
        "strata",
        capture = captures
    )
    model <- simulate(
        data_table,
        "strata",
        suppress,
        method,
        formula_selection,
        opts_stepwise,
        capture = captures
    )

    out <- compare(model, ground_truth)
    return(out)
}
