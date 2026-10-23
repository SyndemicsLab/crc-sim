################################################################################
# File: data_mutating.R                                                        #
# Project: crcsim                                                              #
# Created Date: 2026-10-21                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-21                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' A function used to convert all integer columns to numeric. This is necessary
#' as the estimation functions require numeric data columns.
#'
#' @param df a data frame
#' @returns a data frame with all integer columns converted to numeric
#'
#' @importFrom dplyr mutate across where
#' @keywords internal
#' @export
all_int_cols_to_numeric <- function(df) {
    return(mutate(df, across(where(is.integer), as.numeric)))
}
