################################################################################
# File: test_doubly_robust.R                                                   #
# Project: crcsim                                                              #
# Created Date: 2026-08-31                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("doubly_robust_estimator combines plugin and correction terms", {
    q_1 <- c(0.5, 0.4)
    q_2 <- c(0.5, 0.8)
    q_12 <- c(0.1, 0.2)
    phi <- c(-0.5, 0.5)

    plugin_value <- plugin_estimator(q_1, q_2, q_12)
    result <- crcsim:::doubly_robust_estimator(q_1, q_2, q_12, phi)

    expect_type(result, "double")
    expect_length(result, 1)
    expect_equal(result, plugin_value + mean(phi))
})

test_that("doubly_robust_estimator omits missing correction values", {
    q_1 <- c(0.5, 0.4)
    q_2 <- c(0.5, 0.8)
    q_12 <- c(0.1, 0.2)
    phi <- c(NA_real_, 0.25)

    result <- crcsim:::doubly_robust_estimator(q_1, q_2, q_12, phi)

    expect_equal(
        result,
        plugin_estimator(q_1, q_2, q_12) + mean(phi, na.rm = TRUE)
    )
})

test_that("doubly_robust_estimator enforces a minimum estimate of one", {
    result <- crcsim:::doubly_robust_estimator(
        q_1 = 0.5,
        q_2 = 0.5,
        q_12 = 0.125,
        phi = -10
    )

    expect_equal(result, 1)
})
