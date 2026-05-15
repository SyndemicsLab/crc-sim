################################################################################
# File: crc_plugin.R                                                           #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' CRC with Plugin Estimators
#' @description This function runs CRC using a plugin estimator for the
#' nuisance functions.
#' The plugin estimator is specified in the \code{PluginOptions} object passed
#' to the \code{opts} argument. Supported plugin methods are "logit", "ranger",
#' "logitranger", and "gam". The function will call the appropriate estimation
#' function based on the specified plugin method and return the TMLE estimate
#' of population size and confidence intervals.
#'
#' @param data a data frame containing the capture data. All columns must be
#' numeric type as long as we use the \code{drpop} package.
#' @param opts a \code{PluginOptions} object specifying the plugin method and
#'  parameters for the TMLE estimation.
#' @returns a list containing the TMLE estimate of population size and
#' confidence intervals.
#' @export
crc_plugin <- function(data, opts) {
    if (!inherits(opts, "PluginOptions")) {
        stop(paste(
            "Invalid PluginOptions object provided.",
            "Cannot run CRC with plugin estimator."
        ))
    }
    # TODO: more validation on the data
    data <- all_int_cols_to_numeric(data)
    result <- switch(
        opts$plugin_method,
        logit = logit_estimate(
            data,
            nfolds = opts$nfolds,
            margin = opts$threshold
        ),
        ranger = ranger_estimate(
            data,
            nfolds = opts$nfolds,
            margin = opts$threshold
        ),
        logitranger = logitranger_estimate(
            data,
            nfolds = opts$nfolds,
            margin = opts$threshold
        ),
        gam = gam_estimate(
            data,
            nfolds = opts$nfolds,
            margin = opts$threshold
        ),
        stop(paste(
            "Invalid plugin method specified in PluginOptions object.",
            "Supported methods are 'logit', 'ranger', 'logitranger', and 'gam'."
        ))
    )
    # TODO: check the options estimator and return the applicable result
    # Map should be 1: Doubly Robust, 2: Simple Plugin, 3: TMLE
    return(result)
}
