################################################################################
# File: extract.R                                                              #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-02-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Function to extract the ground truth from data created through the
#' \code{create_data} method
#'
#' @param data_table data.table from \code{create_data}
#' @param capture list: strings of captures
#' @param group character of strata column to extract on
#'
#' @import data.table
#' @returns list the same length of \code{group}
#'
#' @keywords internal
extract_ground_truth <- function(
    data_table,
    group,
    capture = c("APCD", "BSAS", "Casemix", "Death", "Matris", "PMP")
) {
    tmp <- NULL
    data_table <- data_table[, tmp := rowSums(.SD), .SDcols = capture]
    data_table <- data_table[tmp == 0, ][, tmp := NULL]

    out <- as.list(data_table[["N_ID"]])
    names(out) <- ifelse(missing(group), "base", as.list(data_table[[group]]))

    return(out)
}
