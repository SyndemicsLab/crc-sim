################################################################################
# File: test_plugin.R                                                          #
# Project: crcsim                                                              #
# Created Date: 2026-08-31                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("plugin_estimator averages inverse conditional capture", {
    q_1 <- c(0.5, 0.4, 0.8)
    q_2 <- c(0.5, 0.8, 0.5)
    q_12 <- c(0.1, 0.2, 0.4)

    result <- plugin_estimator(
        q_1 = q_1,
        q_2 = q_2,
        q_12 = q_12
    )

    expect_type(result, "double")
    expect_length(result, 1)
    expect_equal(result, mean(1 / conditional_capture(q_1, q_2, q_12)))
})

test_that("plugin_estimator propagates missing conditional probabilities", {
    result <- plugin_estimator(
        q_1 = c(0.5, 0.4),
        q_2 = c(0.5, 0.8),
        q_12 = c(NA_real_, 0.2)
    )

    expect_true(is.na(result))
})

test_that("plugin_estimator accepts scalar nuisance values", {
    result <- plugin_estimator(0.5, 0.4, 0.1)

    expect_equal(result, 2)
})
