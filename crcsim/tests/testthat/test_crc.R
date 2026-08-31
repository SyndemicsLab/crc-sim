################################################################################
# File: test_crc.R                                                             #
# Project: crcsim                                                              #
# Created Date: 2026-08-31                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

make_crc_frequency_fixture <- function() {
    set.seed(123)
    model_data <- create_data(
        n_individuals = 300,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    return(model_data)
}

make_crc_aic_options <- function() {
    formulas <- list(
        stats::as.formula("N_ID ~ capture_1 + capture_2"),
        stats::as.formula("N_ID ~ capture_1 * capture_2")
    )
    aic_options <- AICOptions$new(
        model = "poisson",
        capture_columns = c("capture_1", "capture_2", "capture_3"),
        formula = formulas,
        frequency_col_name = "N_ID"
    )
    return(aic_options)
}

make_crc_stepwise_options <- function() {
    stepwise_options <- StepwiseOptions$new(
        model = "poisson",
        capture_columns = c("capture_1", "capture_2", "capture_3"),
        threshold = 0.05,
        direction = "both",
        frequency_col_name = "N_ID",
        interaction_limit = 2
    )
    return(stepwise_options)
}

make_crc_estimator_options <- function() {
    estimator_options <- EstimatorOptions$new(
        method = "plugin",
        capture_columns = c("capture_1", "capture_2"),
        threshold = 0.01,
        nuisance_function = "logit",
        nfolds = 1
    )
    return(estimator_options)
}

test_that("crc dispatches AIC options to aic_selection", {
    model_data <- make_crc_frequency_fixture()
    opts <- make_crc_aic_options()

    expect_equal(
        quiet_glm_call(crc(model_data, opts)),
        quiet_glm_call(aic_selection(model_data, opts))
    )
})

test_that("crc dispatches stepwise options to stepwise_selection", {
    model_data <- make_crc_frequency_fixture()
    opts <- make_crc_stepwise_options()

    expect_equal(
        quiet_glm_call(crc(model_data, opts)),
        quiet_glm_call(stepwise_selection(model_data, opts))
    )
})

test_that("crc dispatches estimator options to row_level_estimation", {
    model_data <- data.frame(
        capture_1 = c(1, 0, 1, 0, 0),
        capture_2 = c(0, 1, 1, 0, 0)
    )
    opts <- make_crc_estimator_options()

    expect_equal(
        quiet_glm_call(crc(model_data, opts)),
        quiet_glm_call(row_level_estimation(model_data, opts))
    )
})

test_that("crc rejects a bare FrequencyOptions object", {
    opts <- FrequencyOptions$new(
        model = "poisson",
        capture_columns = c("capture_1", "capture_2"),
        threshold = 0.05,
        formulas = NULL,
        frequency_col_name = "N_ID"
    )

    expect_error(
        crc(data.frame(), opts),
        "Invalid FrequencyOptions object provided"
    )
})

test_that("crc rejects objects outside the options hierarchy", {
    expect_error(
        crc(data.frame(), list()),
        "Invalid options object provided"
    )
})
