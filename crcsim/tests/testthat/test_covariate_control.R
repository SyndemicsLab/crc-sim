################################################################################
# File: test_covariate_control_unit.R                                          #
# Project: crcsim                                                              #
# Created Date: 2026-08-28                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("build_covariate_data returns the requested shape and names", {
    result <- build_covariate_data(4, 3)

    expect_s3_class(result, "data.frame")
    expect_equal(dim(result), c(4, 3))
    expect_named(
        result,
        c("covariate_1", "covariate_2", "covariate_3")
    )
})

test_that("build_covariate_data returns binary integer values", {
    set.seed(1)
    result <- build_covariate_data(100, 2)

    expect_type(result$covariate_1, "integer")
    expect_type(result$covariate_2, "integer")
    expect_true(all(unlist(result) %in% c(0L, 1L)))
})

test_that("build_covariate_data is reproducible with the same seed", {
    set.seed(42)
    first <- build_covariate_data(5, 2)
    set.seed(42)
    second <- build_covariate_data(5, 2)

    expect_identical(first, second)
})

test_that("build_covariate_data supports zero covariates", {
    result <- build_covariate_data(4, 0)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 4)
    expect_equal(ncol(result), 0)
    expect_length(names(result), 0)
})

test_that("relocate_covariates moves selected columns to the end", {
    data <- data.frame(
        capture_1 = 1:2,
        covariate_1 = 3:4,
        capture_2 = 5:6,
        covariate_2 = 7:8
    )

    result <- relocate_covariates(data, c("covariate_1", "covariate_2"))

    expect_named(
        result,
        c("capture_1", "capture_2", "covariate_1", "covariate_2")
    )
    expect_equal(result$covariate_1, data$covariate_1)
    expect_equal(result$covariate_2, data$covariate_2)
    expect_equal(result$capture_1, data$capture_1)
    expect_equal(result$capture_2, data$capture_2)
})

test_that("relocate_covariates preserves the requested relative order", {
    data <- data.frame(
        first = 1:2,
        covariate_a = 3:4,
        middle = 5:6,
        covariate_b = 7:8,
        last = 9:10
    )

    result <- relocate_covariates(data, c("covariate_b", "covariate_a"))

    expect_named(
        result,
        c("first", "middle", "last", "covariate_b", "covariate_a")
    )
})

test_that("relocate_covariates returns unchanged data for no columns", {
    data <- data.frame(first = 1:2, second = 3:4)

    expect_identical(relocate_covariates(data, character()), data)
})

test_that("relocate_covariates reports missing columns", {
    data <- data.frame(first = 1:2)

    expect_error(
        relocate_covariates(data, "missing"),
        "Can't select columns that don't exist"
    )
})
