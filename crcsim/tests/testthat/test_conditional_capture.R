################################################################################
# File: test_conditional_capture.R                                             #
# Project: crcsim                                                              #
# Created Date: 2026-08-31                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("conditional_capture returns a good result", {
    q_1 <- 0.4
    q_2 <- 0.5
    q_12 <- 0.2

    result <- crcsim:::conditional_capture(q_1, q_2, q_12)

    expect_type(result, "double")
    expect_equal(result, q_12 / (q_1 * q_2))
})

test_that("conditional_capture handles q_1 being zero", {
    q_1 <- 0
    q_2 <- 0.5
    q_12 <- 0.2

    expect_warning(
        result <- crcsim:::conditional_capture(q_1, q_2, q_12),
        "q_1 or q_2 is zero"
    )

    expect_true(is.na(result))
})

test_that("conditional_capture handles q_2 being zero", {
    q_1 <- 0.4
    q_2 <- 0
    q_12 <- 0.2

    expect_warning(
        result <- crcsim:::conditional_capture(q_1, q_2, q_12),
        "q_1 or q_2 is zero"
    )

    expect_true(is.na(result))
})

test_that("conditional_capture handles both q_1 and q_2 being zero", {
    q_1 <- 0
    q_2 <- 0
    q_12 <- 0.2

    expect_warning(
        result <- crcsim:::conditional_capture(q_1, q_2, q_12),
        "q_1 or q_2 is zero"
    )

    expect_true(is.na(result))
})
