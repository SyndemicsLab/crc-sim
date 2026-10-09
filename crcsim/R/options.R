################################################################################
# File: options.R                                                              #
# Project: crcsim                                                              #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-09-08                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

## NOTE: This file turns off formatting for many of the R6 class definitions to
## preserve the linting of class names. Please pay attention to the formatting
## when editing this file!

#' CRC Options Classes
#' @description Base options shared by CRC estimation workflows.
#' @export
# fmt: skip
Options <- R6::R6Class( # nolint: object_name_linter
    "Options",
    public = list(
        #' @field model Shared model or estimator method field.
        model = NULL,
        #' @field capture_columns Character vector naming the binary capture
        #' indicator columns.
        capture_columns = NULL,
        #' @field threshold Shared threshold or probability margin.
        threshold = NULL,

        #' @description Create a new \code{Options} instance.
        initialize = function(model, capture_columns, threshold) {
            self$model <- model
            self$capture_columns <- capture_columns
            self$threshold <- threshold
            return(self)
        }
    )
)

#' Log-linear Options Class
#' @description This class contains the complete configuration for selecting
#' and fitting log-linear capture-recapture models.
#'
#' @param capture_columns Character vector naming the binary capture
#' indicator columns.
#' @param frequency_col_name Character scalar naming the frequency column.
#' @param model_family Character scalar identifying the log-linear model family,
#' either "poisson" or "negbin".
#' @param selection_method Character scalar identifying the formula selection
#' method, either "aic" or "stepwise".
#' @param selection_criterion Character scalar identifying the information
#' criterion, either "AIC" or "BIC".
#' @param selection_options Named list containing method-specific settings.
#' @export
# fmt: skip
LoglinearOptions <- R6::R6Class( # nolint: object_name_linter
    "LoglinearOptions",
    inherit = Options,
    public = list(
        #' @field frequency_col_name Character scalar naming the frequency
        #' column in the aggregated CRC data.
        frequency_col_name = NULL,
        #' @field model_family Character scalar identifying the log-linear
        #' model family.
        model_family = NULL,
        #' @field selection_method Character scalar identifying the formula
        #' selection method.
        selection_method = NULL,
        #' @field selection_criterion Character scalar identifying the
        #' information criterion.
        selection_criterion = NULL,
        #' @field selection_options Named list of method-specific settings.
        selection_options = NULL,

        #' @description Create a new \code{LoglinearOptions} instance.
        #' @return The initialized \code{LoglinearOptions} object.
        initialize = function(
            capture_columns,
            frequency_col_name = "N_ID",
            model_family = "poisson",
            selection_method = "aic",
            selection_criterion = "AIC",
            selection_options = list()
        ) {
            if (
                !is.null(capture_columns) &&
                    (!is.character(capture_columns) ||
                        anyNA(capture_columns) ||
                        anyDuplicated(capture_columns) > 0)
            ) {
                stop("capture_columns must be unique character names or NULL.")
            }
            if (
                length(frequency_col_name) != 1 ||
                    !is.character(frequency_col_name) ||
                    is.na(frequency_col_name) ||
                    !nzchar(frequency_col_name)
            ) {
                stop("frequency_col_name must be a non-empty character scalar.")
            }
            model_family <- match.arg(model_family, c("poisson", "negbin"))
            selection_method <- match.arg(
                selection_method,
                c("aic", "stepwise")
            )
            selection_criterion <- match.arg(
                selection_criterion,
                c("AIC", "BIC")
            )
            if (
                !is.list(selection_options) ||
                    (
                        length(selection_options) > 0 &&
                            (
                                is.null(names(selection_options)) ||
                                    any(names(selection_options) == "") ||
                                    anyDuplicated(names(selection_options)) > 0
                            )
                    )
            ) {
                stop("selection_options must be a named list.")
            }

            allowed_options <- if (selection_method == "aic") {
                "formulas"
            } else {
                c("direction", "interaction_limit")
            }
            unknown_options <- setdiff(
                names(selection_options),
                allowed_options
            )
            if (length(unknown_options) > 0) {
                stop(
                    "Unsupported selection option(s): ",
                    paste(unknown_options, collapse = ", ")
                )
            }

            if (selection_method == "aic") {
                formulas <- selection_options[["formulas"]]
                if (is.null(formulas)) {
                    if (is.null(capture_columns)) {
                        stop(
                            paste(
                                "capture_columns is required when formulas",
                                "are not provided."
                            )
                        )
                    }
                    formulas <- formula_list(
                        frequency_col_name,
                        capture_columns
                    )
                    selection_options[["formulas"]] <- formulas
                }
            } else {
                direction <- selection_options[["direction"]]
                if (is.null(direction)) {
                    direction <- "both"
                }
                interaction_limit <- selection_options[["interaction_limit"]]
                if (is.null(interaction_limit)) {
                    interaction_limit <- 2
                }
                selection_options[["direction"]] <- match.arg(
                    direction,
                    c("both", "backward", "forward")
                )
                if (
                    length(interaction_limit) != 1 ||
                        !is.numeric(interaction_limit) ||
                        !is.finite(interaction_limit) ||
                        interaction_limit < 1 ||
                        interaction_limit %% 1 != 0
                ) {
                    stop("interaction_limit must be a positive integer.")
                }
                selection_options[["interaction_limit"]] <- as.integer(
                    interaction_limit
                )
            }

            super$initialize(
                model = NULL,
                capture_columns = capture_columns,
                threshold = NULL
            )
            self$frequency_col_name <- frequency_col_name
            self$model_family <- model_family
            self$selection_method <- selection_method
            self$selection_criterion <- selection_criterion
            self$selection_options <- selection_options
            return(self)
        }
    )
)

#' Plugin Options Class
#' @description This class defines the options for the plugin estimator.
#'
#' @param method Character scalar naming the estimation method.
#' @param capture_columns Character vector naming the binary capture indicator
#' columns.
#' @param threshold Numeric scalar giving the probability margin.
#' @param nuisance_function Character scalar identifying the nuisance model.
#' @param nfolds Integer scalar giving the number of cross-validation folds.
#' @export
# fmt: skip
EstimatorOptions <- R6::R6Class( # nolint: object_name_linter
    "EstimatorOptions",
    inherit = Options,
    public = list(
        #' @field nuisance_function Character scalar identifying the function
        #' to use to estimate nuisance parameters.
        nuisance_function = NULL,
        #' @field nfolds Integer scalar giving the number of cross-validation
        #' folds.
        nfolds = NULL,

        #' @description Create a new \code{EstimatorOptions} instance.
        initialize = function(
            method,
            capture_columns,
            threshold,
            nuisance_function,
            nfolds
        ) {
            super$initialize(method, capture_columns, threshold)
            self$nuisance_function <- nuisance_function
            self$nfolds <- nfolds
            return(self)
        }
    )
)
