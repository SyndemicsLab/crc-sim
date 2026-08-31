################################################################################
# File: loglinear.R                                                            #
# Project: crc-sim                                                             #
# Created Date: 2026-05-14                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-06-16                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Fit a log-linear model to the provided data using the specified formula.
#' Generally used as a helper function within the CRC function to fit models
#' for different formula specifications. The user can specify the model family
#' to use, either "poisson" or "negbin".
#'
#' @param data a data frame containing the observed capture histories and
#' frequency column.
#' @param formula_object a formula object specifying the log-linear model to
#' fit.
#' @param model_family a string specifying the model family to use, either
#' "poisson" or "negbin". Default is "poisson".
#' @return a fitted model object resulting from the log-linear model fit.
#'
#' @importFrom stats glm
#' @importFrom MASS glm.nb
#'
#' @export
fit_loglinear_model <- function(
    model_data,
    formula_object,
    model_family = c("poisson", "negbin")
) {
    model_family <- match.arg(model_family)
    if (model_family == "poisson") {
        return(glm(formula_object, data = model_data, family = "poisson"))
    }
    return(glm.nb(formula = formula_object, data = model_data))
}
