################################################################################
# File: crc.R                                                                  #
# Project: crcsim                                                              #
# Created Date: 2026-02-23                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-09-04                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Spatial Capture Re-Capture
#'
#' @description A method for estimation of 'unknowns' through knowledge about
#' knowns as described by Barocas, Joshua A et al. "Estimated Prevalence of
#' Opioid Use Disorder in Massachusetts, 2011-2015: A Capture-Recapture
#' Analysis." doi:10.2105/AJPH.2018.304673
#'
#' This implementation supports poisson and negative binomial regression models
#' with either AIC-based/stepwise formula selection or plugin estimators.
#'
#' @param model_data Dataframe: A dataframe containing a frequency column and
#' binary
#' columns indicating involvement in the given database
#' @param opts Options Object: An object containing the options for the CRC
#' estimation. This should be an instance of one of the options classes defined
#' in options.R, such as LoglinearOptions or EstimatorOptions.
#'
#' @export
crc <- function(model_data, opts, seed = NULL, diagnostics = "quiet") {
    if (inherits(opts, "LoglinearOptions")) {
        if (opts$selection_method == "aic") {
            return(aic_selection(model_data, opts))
        }
        if (opts$selection_method == "stepwise") {
            return(stepwise_selection(model_data, opts))
        }
        stop("Invalid LoglinearOptions selection method.")
    } else if (inherits(opts, "EstimatorOptions")) {
        return(row_level_estimation(model_data, opts, seed, diagnostics))
    }
    stop("Invalid options object provided.")
}
