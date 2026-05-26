################################################################################
# File: nuisance_fit.R                                                         #
# Project: crc-sim                                                             #
# Created Date: 2026-05-20                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-22                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Empirical mean for each capture list and their product for the provided j
#' and k indices. This is used as a fallback when there are no covariates (i.e.
#' when the number of columns in train is less than or equal to n_capture_cols).
#'
#' @param train Training data frame
#' @param test Test data frame
#' @param j Index of first capture list
#' @param k Index of second capture list
#' @return List with elements q1, q2, q12 (predicted capture probabilities)
#'
#' @keywords internal
#' @export
empirical_qhats <- function(train, test, j, k) {
    if (ncol(train) <= max(j, k)) {
        stop(
            paste0(
                "j and k indices must be less than or equal to the number of",
                "columns in train"
            )
        )
    }
    q_j <- mean(train[[j]])
    q_k <- mean(train[[k]])
    q_jk <- mean(train[[j]] * train[[k]])
    return(list(
        q_j = rep(q_j, nrow(test)),
        q_k = rep(q_k, nrow(test)),
        q_jk = rep(q_jk, nrow(test)),
        y_j = test[[j]],
        y_k = test[[k]]
    ))
}

#' Fit a logistic regression model to the capture list given the covariates.
#' This is a helper function for logit_qhats to avoid code repetition.
#'
#' @param data Data frame containing the capture columns and covariates.
#' @param i Index of the capture column to fit.
#' @param n Number of capture columns (the first n columns are capture columns).
#' @return Fitted logistic regression model.
#'
#' @importFrom stats glm binomial
#' @keywords internal
fit_glm <- function(data, formula) {
    return(stats::glm(
        formula,
        data = data,
        family = stats::binomial(link = "logit")
    ))
}

#' Fit a generalized additive model to the capture list given the covariates.
#' This is a helper function for gam_qhats to avoid code repetition.
#'
#' @param data Data frame containing the capture columns and covariates.
#' @param i Index of the capture column to fit.
#' @param n Number of capture columns (the first n columns are capture columns).
#' @return Fitted GAM model.
#'
#' @importFrom mgcv gam
#' @importFrom stats binomial
#' @keywords internal
fit_gam <- function(data, formula) {
    return(mgcv::gam(
        formula,
        data = data,
        family = stats::binomial(link = "logit")
    ))
}

#' Fit a random forest model to the capture list given the covariates.
#' This is a helper function for random_forest_qhats to avoid code repetition.
#'
#' @param data Data frame containing the capture columns and covariates.
#' @param i Index of the capture column to fit.
#' @param n Number of capture columns (the first n columns are capture columns).
#' @param ... Additional arguments passed to ranger().
#' @return Fitted random forest model.
#'
#' @importFrom ranger ranger
#' @keywords internal
fit_random_forest <- function(data, formula) {
    return(ranger::ranger(
        formula,
        data = data,
        probability = TRUE,
        classification = TRUE
    ))
}

#' Predict capture probabilities from a fitted model.
#' This is a helper function for the qhat_* functions to avoid code repetition.
#'
#' @param model Fitted model.
#' @param data Data frame.
#' @returns Predicted probabilities from the specified model.
#'
#' @importFrom stats predict
#' @keywords internal
predict_stats_model_qhat <- function(model, data) {
    return(stats::predict(model, newdata = data, type = "response"))
}

#' Predict capture probabilities from a fitted ranger model.
#' This is a helper function for random_forest_qhats to avoid code repetition.
#' Annoyingly, ranger expects a `data` argument rather than `newdata`, and the
#' predict function for ranger models returns a list with "predictions" element
#' that contains the predicted probabilities, so we need a separate function to
#' extract those probabilities.
#'
#' @param model Fitted ranger model.
#' @param data Data frame.
#' @returns Predicted probabilities from the specified ranger model.
#'
#' @importFrom stats predict
#' @keywords internal
predict_ranger_model_qhat <- function(model, data) {
    return(stats::predict(model, data)$predictions[, "1"])
}

#' Logistic regression nuisance function for capture-recapture
#' Fits logistic regression models on train and predicts capture
#' probabilities on test for a specific list pair (j, k).
#'
#' @param train Training data frame with n_capture_cols capture columns and
#' covariates.
#' @param test Test data frame with n_capture_cols capture columns and
#' covariates.
#' @param n_capture_cols Number of capture lists (columns).
#' @param j Index of first capture list.
#' @param k Index of second capture list.
#' @param margin Margin parameter (not used in logit but kept for interface).
#' @param ... Additional arguments (unused).
#' @returns List with elements q1, q2, q12 (predicted capture probabilities).
#'
#' @importFrom stats formula
#' @keywords internal
#' @export
glm_qhats <- function(
    train,
    test,
    n_capture_cols,
    j,
    k,
    margin = 0.005,
    ...
) {
    # no covariates
    if (ncol(train) <= n_capture_cols) {
        return(empirical_qhats(train, test, j, k))
    }

    # Fit logit models for each list
    j_colname <- colnames(train)[j]
    k_colname <- colnames(train)[k]

    j_form <- stats::formula(paste0(j_colname, "*(1 - ", k_colname, ") ~ ."))
    k_form <- stats::formula(paste0(k_colname, "*(1 - ", j_colname, ") ~ ."))
    jk_form <- stats::formula(paste0(j_colname, "*", k_colname, " ~ ."))

    training_data <- train[,
        c(j, k, (n_capture_cols + 1):ncol(train)),
        drop = FALSE
    ]

    q_j <- fit_glm(training_data, j_form) |>
        predict_stats_model_qhat(test) |>
        pmax(margin)

    q_k <- fit_glm(training_data, k_form) |>
        predict_stats_model_qhat(test) |>
        pmax(margin)

    q_jk <- fit_glm(training_data, jk_form) |>
        predict_stats_model_qhat(test) |>
        pmax(margin)

    return(list(
        q_j = q_j,
        q_k = q_k,
        q_jk = q_jk,
        y_j = test[[j]],
        y_k = test[[k]]
    ))
}

#' GAM nuisance function for capture-recapture
#' Fits generalized additive models on train and predicts capture
#' probabilities on test for a specific list pair (j, k).
#'
#' @param train Training data frame with n_capture_cols capture columns and
#' covariates.
#' @param test Test data frame with n_capture_cols capture columns and
#' covariates.
#' @param n_capture_cols Number of capture lists (columns).
#' @param j Index of first capture list.
#' @param k Index of second capture list.
#' @param margin Margin parameter (not used in GAM but kept for interface).
#' @param ... Additional arguments (unused).
#' @returns List with elements q1, q2, q12 (predicted capture probabilities).
#'
#' @importFrom stats formula
#' @keywords internal
#' @export
gam_qhats <- function(train, test, n_capture_cols, j, k, margin = 0.005, ...) {
    if (ncol(train) <= n_capture_cols) {
        return(empirical_qhats(train, test, j, k))
    }

    # Fit logit models for each list
    j_colname <- colnames(train)[j]
    k_colname <- colnames(train)[k]
    cov_colnames <- colnames(train)[-(1:n_capture_cols)]

    j_form <- stats::formula(paste0(
        j_colname,
        " ~ ",
        paste(cov_colnames, collapse = " + ")
    ))
    k_form <- stats::formula(paste0(
        k_colname,
        " ~ ",
        paste(cov_colnames, collapse = " + ")
    ))
    jk_form <- stats::formula(paste0(
        j_colname,
        "*",
        k_colname,
        " ~ ",
        paste(cov_colnames, collapse = " + ")
    ))

    # Fit GAM models for each list
    q_j <- fit_gam(train, j_form) |>
        predict_stats_model_qhat(test) |>
        pmax(margin)
    q_k <- fit_gam(train, k_form) |>
        predict_stats_model_qhat(test) |>
        pmax(margin)
    q_jk <- fit_gam(train, jk_form) |>
        predict_stats_model_qhat(test) |>
        pmax(margin)

    return(list(
        q_j = q_j,
        q_k = q_k,
        q_jk = q_jk,
        y_j = test[[j]],
        y_k = test[[k]]
    ))
}


#' Ranger nuisance function for capture-recapture
#' Fits random forest models (ranger) on train and predicts capture
#' probabilities on test for a specific list pair (j, k).
#'
#' @param train Training data frame with n_capture_cols capture columns and
#' covariates.
#' @param test Test data frame with n_capture_cols capture columns and
#' covariates.
#' @param n_capture_cols Number of capture lists (columns).
#' @param j Index of first capture list.
#' @param k Index of second capture list.
#' @param margin Margin parameter (not used in ranger but kept for interface).
#' @param ... Additional arguments passed to ranger().
#' @returns List with elements q1, q2, q12 (predicted capture probabilities).
#'
#' @importFrom stats formula
#' @keywords internal
#' @export
random_forest_qhats <- function(
    train,
    test,
    n_capture_cols,
    j,
    k,
    margin = 0.005,
    ...
) {
    if (ncol(train) <= n_capture_cols) {
        return(empirical_qhats(train, test, j, k))
    }

    # Fit logit models for each list
    j_colname <- colnames(train)[j]
    k_colname <- colnames(train)[k]

    j_form <- stats::formula(paste0("factor(", j_colname, ") ~ ."))
    k_form <- stats::formula(paste0("factor(", k_colname, ") ~ ."))
    jk_form <- stats::formula(paste0(
        "factor(",
        j_colname,
        "*",
        k_colname,
        ") ~ ."
    ))

    cc <- 1:n_capture_cols
    j_data <- train[, -cc[-j]]
    k_data <- train[, -cc[-k]]
    jk_data <- train[, -cc[-c(j, k)]]

    q_j <- fit_random_forest(j_data, j_form) |>
        predict_ranger_model_qhat(test) |>
        pmax(margin)
    q_k <- fit_random_forest(k_data, k_form) |>
        predict_ranger_model_qhat(test) |>
        pmax(margin)
    q_jk <- fit_random_forest(jk_data, jk_form) |>
        predict_ranger_model_qhat(test) |>
        pmax(margin)

    return(list(
        q_j = q_j,
        q_k = q_k,
        q_jk = q_jk,
        y_j = test[[j]],
        y_k = test[[k]]
    ))
}

# TODO: This can be expanded to include ensemble methods like Super Learner and Ranger-Logit, but we should fix the estimators before adding more nuisance functions to avoid scope creep.
