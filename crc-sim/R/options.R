################################################################################
# File: options.R                                                              #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Your Company                                              #
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
        #' @field model Character scalar identifying either the log-linear
        #' model or the nuisance function to utilize.
        model = NULL,
        #' @field threshold Numeric scalar giving the threshold applied by
        #' threshold-based selection methods or the margin desired by stepwise
        #' selection.
        threshold = NULL,

        #' @description Create a new \\code{Options} instance.
        #' @param model Character scalar identifying either the log-linear
        #' model or the nuisance function to utilize.
        #' @param threshold Numeric scalar giving the threshold applied by
        #' threshold-based selection methods or the margin desired by stepwise
        #' selection.
        #' @return The initialized \\code{Options} object.
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
#' @param formulas Optional formula object or list of formulas to initialize
#' the object with.
#' @export
# fmt:skip
FrequencyOptions <- R6::R6Class( # nolint: object_name_linter
    "FrequencyOptions",
    inherit = Options,
    public = list(
        #' @field formulas Formula object or collection of formulas used when
        #' fitting frequency-based CRC models.
        formulas = NULL,
        #' @field frequency_col_name Character scalar naming the frequency
        #' column.
        frequency_col_name = NULL,

        #' @description Create a new \\code{FrequencyOptions} instance.
        #' @param model Character scalar identifying the model family.
        #' @param threshold Numeric scalar giving the threshold applied by
        #' threshold-based selection methods.
        #' @param formulas Optional formula object or list of formulas used to
        #' initialize the object.
        #' @param frequency_col_name Character scalar naming the frequency
        #' column in the aggregated CRC data.
        #' @return The initialized \\code{FrequencyOptions} object.
        initialize = function(
            model,
            threshold,
            formulas,
            frequency_col_name
        ) {
            super$initialize(model, threshold)
            self$formulas <- formulas
            self$frequency_col_name <- frequency_col_name
            return(self)
        },

        #' @description Append a single formula to \\code{self$formulas}.
        #' @param formula A formula object to append.
        #' @return The updated \\code{FrequencyOptions} object.
        add_formula = function(formula) {
            if (is.null(self$formulas)) {
                self$formulas <- formula
            } else {
                self$formulas <- c(self$formulas, formula)
            }
            return(self)
        },

        #' @description Append multiple formulas to \\code{self$formulas}.
        #' @param formula_list A list of formula objects to append.
        #' @return The updated \\code{FrequencyOptions} object.
        add_formulas = function(formula_list) {
            if (is.null(self$formulas)) {
                self$formulas <- formula_list
            } else {
                self$formulas <- c(self$formulas, formula_list)
            }
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
#' @param interaction_limit Integer scalar giving the maximum order of
#' interactions to include in the stepwise search.
#' @export
#fmt: skip
StepwiseOptions <- R6::R6Class( # nolint: object_name_linter
    "StepwiseOptions",
    inherit = FrequencyOptions,
    public = list(
        #' @field direction Character scalar specifying the stepwise search
        #' direction.
        direction = NULL,
        #' @field capture_indicators Character vector naming the capture history
        #' indicator columns.
        capture_indicators = NULL,
        #' @field interaction_limit Integer scalar giving the maximum order of
        #' interactions to include in the search.
        interaction_limit = 2,

        #' @description Create a new \\code{StepwiseOptions} instance.
        #' @param model Character scalar identifying the model family.
        #' @param threshold Numeric scalar giving the stepwise inclusion
        #' threshold.
        #' @param direction Character scalar specifying the stepwise search
        #' direction.
        #' @param frequency_column Character scalar naming the frequency column
        #' in the aggregated CRC data.
        #' @param capture_indicators Character vector naming the capture history
        #' indicator columns.
        #' @param interaction_limit Integer scalar giving the maximum order of
        #' interactions to include in the search.
        #' @return The initialized \\code{StepwiseOptions} object.
        initialize = function(
            model,
            threshold,
            direction,
            frequency_col_name = "N_ID",
            capture_indicators = NULL,
            interaction_limit = 2
        ) {
            super$initialize(
                model,
                threshold,
                formulas = NULL,
                frequency_col_name = frequency_col_name
            )
            self$direction <- direction
            self$capture_indicators <- capture_indicators
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
    inherit = FrequencyOptions,
    public = list(
        #' @field formulas Formula object or collection of formulas evaluated by
        #' the AIC-based selection routine.
        formulas = NULL,
        #' @description Create a new \\code{AICOptions} instance.
        #' @param model Character scalar identifying the model family.
        #' @param formula Optional formula object to evaluate directly.
        #' @param frequency_col_name Character scalar naming the frequency
        #' column
        #' in the aggregated CRC data.
        #' @param capture_indicators Character vector naming the capture history
        #' indicator columns.
        #' @return The initialized \\code{AICOptions} object.
        initialize = function(
            model,
            formula = NULL,
            frequency_col_name = "N_ID",
            capture_indicators = NULL
        ) {
            super$initialize(
                model,
                threshold = NULL,
                formulas = NULL,
                frequency_col_name = frequency_col_name
            )
            self$formulas <- private$validate_formula_input(
                frequency_col_name,
                capture_indicators,
                formula
            )
            return(self)
        }
    ),
    private = list(
        # Note: These are not Roxygen2 comments because it errors on private
        # members of classes being commented
        # @description Validate or generate the formula input used by
        # \\code{AICOptions}.
        # @param freq_column Character scalar naming the frequency column in
        # the aggregated CRC data.
        # @param binary_variables Character vector naming the binary capture
        # indicator columns.
        # @param formula Optional formula object supplied by the caller.
        # @return A formula object or collection of formulas for AIC model
        # evaluation.
        validate_formula_input = function(
            frequency_col_name,
            binary_variables,
            formula
        ) {
            if (!is.null(formula)) {
                return(formula)
            }
            if (is.null(frequency_col_name) || is.null(binary_variables)) {
                stop(
                    paste(
                        "If formula is not provided, frequency_col_name and",
                        "binary_variables must be specified."
                    )
                )
            }
            return(formula_list(frequency_col_name, binary_variables))
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
        #' @field list_pair Character vector of length two naming the capture
        #' lists used by the plugin estimator.
        list_pair = NULL,
        #' @field nfolds Integer scalar giving the number of cross-validation
        #' folds.
        nfolds = NULL,
        #' @field plugin_estimator Character scalar naming the plugin estimator
        #' implementation.
        plugin_estimator = NULL,

        #' @description Create a new \\code{PlugInOptions} instance.
        #' @param model Character scalar identifying the model family.
        #' @param threshold Numeric scalar giving the threshold applied by
        #' threshold-based selection methods.
        #' @param list_pair Character vector of length two naming the capture
        #' lists used by the plugin estimator.
        #' @param nfolds Integer scalar giving the number of cross-validation
        #' folds.
        #' @param plugin_estimator Character scalar naming the plugin estimator
        #' implementation.
        #' @return The initialized \\code{PlugInOptions} object.
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
