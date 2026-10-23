################################################################################
# File: estimator_variants.R                                                   #
# Project: crcsim                                                              #
# Created Date: 2026-10-23                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' @title Doubly Robust Estimator Class
#'
#' @description Creates a doubly robust estimator object with the specified
#' parameters. This estimator object is used for the DR estimator as well as
#' the TMLE estimator.
#'
#' @details This is meant to be used in conjunction with the `drpop` R package
#' for doubly robust estimation. It provides a structured way to ensure all the
#' necessary parameters are applied. Possible estimator functions include:
#'
#' - "logit": Logistic regression
#' - "gam": Generalized additive model
#' - "ranger": Random forest via the `ranger` package
#' - "mlogit": Multinomial logistic regression
#' - "rangerlogit": Combination of random forest and logistic regression
#' - "sl": Super learner
#'
#' @param n_capture_cols Number of capture columns in the dataset.
#' @param estimator_func Estimator function to use. One of "logit", "gam",
#' "ranger", "mlogit", "rangerlogit", or "sl". Default is "logit".
#' @param n_folds Number of folds for cross-validation (default is 5).
#' @param margin Margin parameter for the estimator (default is 0.005).
#' @param diagnostics Controls GLM fit diagnostics: "quiet" (default),
#' "warning", or "verbose".
#'
#' @return A doubly robust estimator object of class "doubly_robust_estimator"
#' and "estimator".
#'
#' @export
doubly_robust_estimator <- function(
    estimator_func = "logit",
    n_folds = 5,
    margin = 0.005,
    diagnostics = "quiet"
) {
    estimator_func <- match.arg(
        estimator_func,
        c("logit", "gam", "ranger", "mlogit", "rangerlogit", "sl")
    )
    diagnostics <- match.arg(diagnostics, c("quiet", "warning", "verbose"))
    estim <- list(
        estimator_func = estimator_func,
        n_folds = n_folds,
        margin = margin,
        diagnostics = diagnostics
    )
    class(estim) <- c("doubly_robust_estimator", "estimator")
    return(estim)
}

#' @title Log Linear Estimator
#'
#' @description A log linear estimator object for modeling count data using
#' specified model family and selection criteria.
#'
#' @details The log linear estimator allows for flexible modeling of count data.
#' The `model_family` parameter specifies the distribution family, either
#' "poisson" or "negbin".
#' The `selection` parameter determines the model selection method, either
#' "aic" or "stepwise".
#' The `criterion` parameter specifies the information criterion to use, either
#' "AIC" or "BIC".
#' The `freq_col_name` parameter indicates the column name for frequency counts.
#'
#' @param n_capture_cols Number of capture columns in the dataset.
#' @param model_family Distribution family for the model (default is "poisson").
#' @param selection Model selection method, either "aic" or "stepwise" (default
#' is "aic").
#' @param criterion Information criterion to use, either "AIC" or "BIC"
#' (default is "AIC").
#' @param freq_col_name Column name for frequency counts (default is "N_ID").
#' @param selection_opts List of options for model selection (default is an
#' empty list). Optional parameters include:
#' - `direction`: Direction of stepwise selection, either "both", "forward", or
#'  "backward" (default is "both").
#' - `interaction_limit`: Maximum number of interactions allowed in the model
#'  (default is 1).
#' - `formulas`: A list of formula objects specifying the models to fit and
#'  compare using AIC.
#'
#' @return A log linear estimator object of class "log_linear_estimator" and
#' "estimator".
#'
#' @export
log_linear_estimator <- function(
    model_family = "poisson",
    selection = "aic",
    criterion = "AIC",
    freq_col_name = "N_ID",
    selection_opts = list()
) {
    model_family <- match.arg(model_family, c("poisson", "negbin"))
    selection <- match.arg(selection, c("aic", "stepwise"))
    criterion <- match.arg(criterion, c("AIC", "BIC"))

    estim <- list(
        n_capture_cols = n_capture_cols,
        model_family = model_family,
        selection = selection,
        criterion = criterion,
        freq_col_name = freq_col_name,
        selection_opts = selection_opts
    )
    class(estim) <- c("log_linear_estimator", "estimator")
    return(estim)
}
