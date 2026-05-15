################################################################################
# File: options.R                                                              #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

## NOTE: This file turns off formatting for many of the R6 class definitions to
## preserve the linting of class names. Please pay attention to the formatting
## when editing this file!

#' CRC Options Classes
#' @description This file defines the R6 classes for the options used in the
#' different CRC formula selection methods. These classes encapsulate the
#' parameters and settings for each method, allowing for organized and
#' structured handling of the options when running CRC with different formula
#' selection approaches.
#'
#' @param model The model family or function to apply to the data. If using a
#' log-linear model, this should be either "poisson" or "negbin". If using
#' a plugin estimator, this should be a "logit", "ranger", "rangerlogit", or
#' "gam".
#' @param threshold The p-value threshold for variable inclusion in the model.
#' @export
# fmt: skip
Options <- R6::R6Class( # nolint: object_name_linter
    "Options",
    public = list(
        model = NULL,
        threshold = NULL,

        initialize = function(model, threshold) {
            self$model <- model
            self$threshold <- threshold
            return(self)
        }
    )
)

#' Frequency Options Class
#' @description This class defines the options for the frequency-based formula
#' selection methods in CRC. It inherits from the base Options class and
#' includes additional parameters specific to the frequency-based methods
#'
#' @param formula A formula object specifying the log-linear model to fit. If
#' NULL, the formula will be determined based on the specified formula
#' selection method.
#' @param frequency_column The name of the column in the data that contains the
#' frequency counts for each capture history.
#' @param capture_indicators A vector of the names of the columns in the data
#' that indicate the capture history (i.e., which lists captured each
#' individual).
#' @export
# fmt:skip
FrequencyOptions <- R6::R6Class( # nolint: object_name_linter
    "FrequencyOptions",
    inherit = Options,
    public = list(
        formulas = NULL,
        frequency_column = NULL,
        capture_indicators = NULL,

        initialize = function(
            model,
            threshold,
            frequency_column = NULL,
            capture_indicators = NULL
        ) {
            super$initialize(model, threshold)
            self$frequency_column <- frequency_column
            self$capture_indicators <- capture_indicators
            return(self)
        }
    )
)

#' Stepwise Options Class
#' @description This class defines the options for the stepwise formula
#' selection method in CRC. It inherits from the base Options class and
#' includes additional parameters specific to the stepwise method.
#'
#' @param direction The direction of the stepwise selection, either "forward",
#' "backward", or "both".
#' @param frequency_column The name of the column in the data that contains the
#' frequency counts for each capture history.
#' @param capture_indicators A vector of the names of the columns in the data
#' that indicate the capture history (i.e., which lists captured each
#' individual).
#' @export
#fmt: skip
StepwiseOptions <- R6::R6Class( # nolint: object_name_linter
    "StepwiseOptions",
    inherit = FrequencyOptions,
    public = list(
        direction = NULL,
        threshold = NULL,
        frequency_column = NULL,
        capture_indicators = NULL,
        interaction_limit = 2,

        initialize = function(
            model,
            threshold,
            direction,
            formula = NULL,
            frequency_column = NULL,
            capture_indicators = NULL,
            interaction_limit = 2
        ) {
            super$initialize(
                model,
                threshold,
                formula,
                frequency_column,
                capture_indicators
            )
            self$direction <- direction
            self$interaction_limit <- interaction_limit
            return(self)
        }
    )
)

#' AIC Options Class
#' @description This class defines the options for the AIC-based formula
#' selection method in CRC. It inherits from the base Options class and
#' includes additional parameters specific to the AIC-based method.
#'
#' @param formula A formula object specifying the log-linear model to fit. If
#' NULL, the formula will be determined based on the specified formula
#' selection method.
#' @param frequency_column The name of the column in the data that contains the
#' frequency counts for each capture history.
#' @param capture_indicators A vector of the names of the columns in the data
#' that indicate the capture history (i.e., which lists captured each
#' individual).
#' @export
# fmt: skip
AICOptions <- R6::R6Class( # nolint: object_name_linter
    "AICOptions",
    inherit = Options,
    public = list(
        formulas = NULL,
        initialize = function(
            model,
            frequency_column = NULL,
            capture_indicators = NULL,
            formula = NULL
        ) {
            super$initialize(
                model,
                threshold = NULL,
                frequency_column,
                capture_indicators,
                formula
            )
            formulas <- private$validate_formula_input(
                frequency_column,
                capture_indicators,
                formula
            )
            return(self)
        }
    ),
    private = list(
        validate_formula_input = function(
            freq_column,
            binary_variables,
            formula
        ) {
            if (!is.null(formula)) {
                return(formula)
            }
            if (is.null(freq_column) || is.null(binary_variables)) {
                stop(
                    paste(
                        "If formula is not provided, frequency_column and",
                        "binary_variables must be specified."
                    )
                )
            }
            return(formula_list(freq_column, binary_variables))
        }
    )
)

#' Plugin Options Class
#' @description This class defines the options for the plugin estimator formula
#' selection method in CRC. It inherits from the base Options class and
#' includes additional parameters specific to the plugin method.
#'
#' @param list_pair A vector of two column names in the data that indicate the
#' two capture lists to be used for the plugin estimator.
#' @param nfolds The number of folds to use for cross-validation when fitting
#' the plugin estimator model.
#' @param plugin_estimator The type of plugin estimator to use, either "logit",
#' "ranger", "rangerlogit", or "gam".
#' @export
# fmt: skip
PlugInOptions <- R6::R6Class( # nolint: object_name_linter
    "PlugInOptions",
    inherit = Options,
    public = list(
        list_pair = NULL,
        nfolds = NULL,
        plugin_estimator = NULL,

        initialize = function(
            model,
            threshold,
            list_pair,
            nfolds,
            plugin_estimator
        ) {
            super$initialize(model, threshold)
            self$list_pair <- list_pair
            self$nfolds <- nfolds
            self$plugin_estimator <- plugin_estimator
            return(self)
        }
    )
)
