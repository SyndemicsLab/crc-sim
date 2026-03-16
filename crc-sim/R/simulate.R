################################################################################
# File: simulate.R                                                             #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-02-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Attempts recovery of the ground truth - fundamentally a wrapper for crc
#'
#' @param data_table data.table from the \code{create.data} step
#' @param capture list: strings of captures
#' @param group stratification string to group by
#' @param suppress numeric: maximum value to suppress at
#' @param method string: selection for spatial capture recapture model - either
#' "poisson" or "negbin"
#' @param formula_selection string: selection for formula decision - either
#' "aic", "corr", or "stepwise"
#' @param opts_stepwise list: list containing 'direction' of 'forward'
#' 'backward' or 'both', and 'threshold': p value threshold for stepwise
#' selection
#'
#' @keywords internal
simulate <- function(
    data_table,
    group,
    suppress,
    method,
    formula_selection,
    opts_stepwise,
    capture = c("APCD", "BSAS", "Casemix", "Death", "Matris", "PMP")
) {
    N_ID <- NULL
    tmp <- NULL

    data_table <- data_table[, tmp := rowSums(.SD), .SDcols = capture][
        tmp != 0,
    ][,
        tmp := NULL
    ]

    if (suppress) {
        data_table <- data_table[!N_ID %in% 1:suppress, ]
    }

    if (!missing(group)) {
        dt_list <- c()
        n_groups <- sort(unique(data_table[[group]]))

        for (i in seq_along(n_groups)) {
            dt_list[[i]] <- data_table[get(group) == n_groups[i], ][,
                paste0(group) := NULL
            ]
        }
        out_list <- lapply(dt_list, function(x) {
            crc(
                x,
                "N_ID",
                binary.variables = capture,
                method,
                formula_selection,
                opts_stepwise = opts_stepwise
            )
        })
        return(out_list)
    } else {
        out <- crc(
            data_table,
            "N_ID",
            binary.variables = capture,
            method,
            formula_selection,
            opts_stepwise = opts_stepwise
        )
        return(out)
    }
}
