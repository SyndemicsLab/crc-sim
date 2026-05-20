################################################################################
# File: test_crc_plugin_logit.R                                                #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-19                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

make_plugin_logit_fixture <- function(n_individuals = 300, n_captures = 3) {
    set.seed(123)
    return(create_data(
        n_individuals = n_individuals,
        n_captures = n_captures,
        return_frequency_table = FALSE
    ))
}

make_plugin_logit_options <- function(estimator = "TMLE") {
    return(EstimatorOptions$new(
        model = "logit",
        threshold = 0.005,
        list_pair = c("capture_1", "capture_2"),
        nfolds = 2,
        estimator = estimator
    ))
}

mock_plugin_result <- function(
    n = c(100, 110, 120),
    cin_l = c(90, 100, 110),
    cin_u = c(110, 120, 130)
) {
    return(list(result = list(n = n, cin.l = cin_l, cin.u = cin_u)))
}

test_that("crc_plugin logit returns estimate and confidence interval", {
    data <- make_plugin_logit_fixture()
    opts <- make_plugin_logit_options(estimator = "TMLE")

    out <- with_mocked_bindings(
        crc_plugin(data, opts),
        logit_estimate = function(data, nfolds, margin) {
            return(mock_plugin_result())
        }
    )

    expect_type(out, "list")
    expect_named(
        out,
        c("estimate", "lower_ci", "upper_ci", "estimator")
    )
    expect_type(out$estimate, "double")
    expect_type(out$lower_ci, "double")
    expect_type(out$upper_ci, "double")
    expect_equal(out$estimator, "TMLE")
    expect_equal(out$estimate, 120)
    expect_equal(out$lower_ci, 110)
    expect_equal(out$upper_ci, 130)
    expect_lte(out$lower_ci, out$estimate)
    expect_lte(out$estimate, out$upper_ci)
})

test_that("crc_plugin logit rejects frequency-table input", {
    set.seed(456)
    data <- create_data(
        n_individuals = 300,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts <- make_plugin_logit_options()

    expect_error(
        crc_plugin(data, opts),
        "Data must not be a frequency table"
    )
})

test_that("crc_plugin logit rejects non-EstimatorOptions opts", {
    data <- make_plugin_logit_fixture()

    expect_error(
        crc_plugin(data, list()),
        "Invalid EstimatorOptions object provided"
    )
})

test_that("crc_plugin logit output has no silent NA branch", {
    data <- make_plugin_logit_fixture()
    opts <- make_plugin_logit_options(estimator = "DR")

    expect_error(
        with_mocked_bindings(
            crc_plugin(data, opts),
            logit_estimate = function(data, nfolds, margin) {
                return(mock_plugin_result(
                    n = c(NA_real_, 110, 120),
                    cin_l = c(90, 100, 110),
                    cin_u = c(110, 120, 130)
                ))
            }
        ),
        "Plugin estimator produced NA output values"
    )
})

test_that("crc_plugin rejects unsupported plugin method", {
    data <- make_plugin_logit_fixture()
    opts <- make_plugin_logit_options()
    opts$model <- "not_a_method"

    expect_error(
        crc_plugin(data, opts),
        "Invalid plugin method specified"
    )
})

test_that("crc_plugin logit propagates estimator errors", {
    data <- make_plugin_logit_fixture()
    opts <- make_plugin_logit_options()

    expect_error(
        with_mocked_bindings(
            crc_plugin(data, opts),
            logit_estimate = function(data, nfolds, margin) {
                stop("synthetic logit failure")
            }
        ),
        "synthetic logit failure"
    )
})
