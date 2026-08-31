################################################################################
# File: test_row_level_estimation.R                                            #
# Project: crcsim                                                              #
# Created Date: 2026-08-31                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

make_estimator_options <- function(
    capture_columns = c("capture_1", "capture_2")
) {
    estimator_options <- EstimatorOptions$new(
        method = "plugin",
        capture_columns = capture_columns,
        threshold = 0.01,
        nusiance_function = "logit",
        nfolds = 1
    )
    return(estimator_options)
}

test_that("row_level_estimation rejects invalid options", {
    data <- data.frame(capture_1 = c(1, 0), capture_2 = c(0, 1))

    expect_error(
        row_level_estimation(data, list()),
        "Invalid EstimatorOptions object provided"
    )
})

test_that("row_level_estimation rejects frequency tables", {
    data <- data.frame(
        capture_1 = c(1, 0),
        capture_2 = c(0, 1),
        N_ID = c(3, 2)
    )

    expect_error(
        row_level_estimation(data, make_estimator_options()),
        "Data must not be a frequency table"
    )
})

test_that("row_level_estimation uses specified capture columns", {
    data <- data.frame(
        capture_1 = c(1, 0, 1, 0, 0),
        capture_2 = c(0, 1, 1, 0, 0),
        binary_covariate = c(0, 1, 0, 1, 1)
    )

    result <- row_level_estimation(
        data,
        make_estimator_options(c("capture_1", "capture_2"))
    )

    expect_named(result, c("listpair", "n", "sigma_n", "ci_l", "ci_u"))
    expect_equal(result$listpair, "1,2")
    expect_equal(result$n, 3)
})

test_that("row_level_estimation infers binary capture columns", {
    data <- data.frame(
        capture_1 = c(1, 0, 1, 0, 0),
        capture_2 = c(0, 1, 1, 0, 0),
        continuous_covariate = c(0.2, 0.4, 0.6, 0.8, 1.0)
    )

    result <- suppressWarnings(row_level_estimation(
        data,
        make_estimator_options(capture_columns = NULL)
    ))

    expect_equal(result$listpair, "1,2")
    expect_equal(result$n, 3)
})

test_that("row_level_estimation removes uncaptured records before estimation", {
    data <- data.frame(
        capture_1 = c(1, 0, 1, 0, 0),
        capture_2 = c(0, 1, 1, 0, 0)
    )

    result <- row_level_estimation(data, make_estimator_options())

    expect_equal(result$listpair, "1,2")
    expect_equal(result$n, 4)
    expect_gte(result$n, result$ci_l)
    expect_gte(result$ci_u, result$n)
})
