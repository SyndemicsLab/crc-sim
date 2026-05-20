################################################################################
# File: test_crc_plugin_ranger.R                                               #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-19                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

make_plugin_ranger_fixture <- function(n_individuals = 300, n_captures = 3) {
    set.seed(123)
    return(create_data(
        n_individuals = n_individuals,
        n_captures = n_captures,
        return_frequency_table = FALSE
    ))
}

make_plugin_ranger_options <- function(estimator = "TMLE") {
    return(EstimatorOptions$new(
        model = "ranger",
        threshold = 0.005,
        list_pair = c("capture_1", "capture_2"),
        nfolds = 2,
        estimator = estimator
    ))
}

mock_plugin_result_ranger <- function(
    n = c(200, 210, 220),
    cin_l = c(190, 200, 210),
    cin_u = c(210, 220, 230)
) {
    return(list(result = list(n = n, cin.l = cin_l, cin.u = cin_u)))
}

test_that("crc_plugin ranger returns estimate and confidence interval", {
    data <- make_plugin_ranger_fixture()
    opts <- make_plugin_ranger_options(estimator = "TMLE")

    out <- with_mocked_bindings(
        crc_plugin(data, opts),
        ranger_estimate = function(data, nfolds, margin) {
            return(mock_plugin_result_ranger())
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
    expect_equal(out$estimate, 220)
    expect_equal(out$lower_ci, 210)
    expect_equal(out$upper_ci, 230)
    expect_lte(out$lower_ci, out$estimate)
    expect_lte(out$estimate, out$upper_ci)
})

test_that("crc_plugin ranger rejects frequency-table input", {
    set.seed(456)
    data <- create_data(
        n_individuals = 300,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts <- make_plugin_ranger_options()

    expect_error(
        crc_plugin(data, opts),
        "Data must not be a frequency table"
    )
})

test_that("crc_plugin ranger output has no silent NA branch", {
    data <- make_plugin_ranger_fixture()
    opts <- make_plugin_ranger_options(estimator = "PI")

    expect_error(
        with_mocked_bindings(
            crc_plugin(data, opts),
            ranger_estimate = function(data, nfolds, margin) {
                return(mock_plugin_result_ranger(
                    n = c(200, NA_real_, 220),
                    cin_l = c(190, 200, 210),
                    cin_u = c(210, 220, 230)
                ))
            }
        ),
        "Plugin estimator produced NA output values"
    )
})

test_that("crc_plugin ranger propagates estimator errors", {
    data <- make_plugin_ranger_fixture()
    opts <- make_plugin_ranger_options()

    expect_error(
        with_mocked_bindings(
            crc_plugin(data, opts),
            ranger_estimate = function(data, nfolds, margin) {
                stop("synthetic ranger failure")
            }
        ),
        "synthetic ranger failure"
    )
})
