################################################################################
# File: test_data_control.R                                                    #
# Project: crc-sim                                                             #
# Created Date: 15 May 2026                                                    #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-06-16                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

# ==============================================================================
# Group 1: create_data() - Row-level dataset
# ==============================================================================

test_that("create_data returns a tibble with n_individuals rows", {
    result <- create_data(n_individuals = 50, n_captures = 2)

    expect_s3_class(result, "tbl_df")
    expect_equal(nrow(result), 50)
})

test_that("create_data capture columns contain only binary values", {
    result <- create_data(n_individuals = 100, n_captures = 3)

    capture_cols <- grep("^capture_", names(result), value = TRUE)
    capture_values <- unlist(result[, capture_cols])

    expect_true(all(capture_values %in% c(0L, 1L)))
})

test_that(
    paste(
        "create_data does not include N_ID when return",
        "frequency_table is FALSE"
    ),
    {
        result <- create_data(
            n_individuals = 50,
            n_captures = 2,
            return_frequency_table = FALSE
        )

        expect_false("N_ID" %in% names(result))
    }
)

# ==============================================================================
# Group 2: create_data() - Aggregated dataset
# ==============================================================================

test_that(
    paste(
        "create_data includes N_ID column when return",
        "frequency_table is TRUE"
    ),
    {
        result <- create_data(
            n_individuals = 100,
            n_captures = 2,
            return_frequency_table = TRUE
        )

        expect_true("N_ID" %in% names(result))
    }
)

test_that("create_data frequency table N_ID sums to n_individuals", {
    result <- create_data(
        n_individuals = 200,
        n_captures = 3,
        return_frequency_table = TRUE
    )

    expect_equal(sum(result$N_ID), 200)
})

test_that("create_data frequency table has all N_ID values >= 1", {
    result <- create_data(
        n_individuals = 200,
        n_captures = 3,
        return_frequency_table = TRUE
    )

    expect_true(all(result$N_ID >= 1))
})

test_that("create_data frequency table has fewer rows than n_individuals", {
    set.seed(2026)
    result <- create_data(
        n_individuals = 500,
        n_captures = 3,
        return_frequency_table = TRUE
    )

    expect_lt(nrow(result), 500)
})

# ==============================================================================
# Group 3: build_individual_records()
# ==============================================================================

test_that("build_individual_records returns correct tibble", {
    covariate_specs <- data.frame(
        distribution = "uniform",
        p1 = 0,
        p2 = 1,
        dtype = "integer"
    )
    result <- build_individual_records(
        n_individuals = 75,
        capture_probs = c(0.5, 0.5),
        covariate_specs = covariate_specs
    )
    expected_names <- c(
        "capture_1",
        "capture_2",
        "covariate_1"
    )

    expect_s3_class(result, "tbl_df")
    expect_equal(nrow(result), 75)
    expect_named(result, expected_names)
})

# ==============================================================================
# Group 4: build_contingency_table()
# ==============================================================================

test_that("build_contingency_table adds an N_ID column", {
    records <- tibble::tibble(
        capture_1 = c(1L, 0L, 1L),
        capture_2 = c(0L, 1L, 1L)
    )
    result <- build_contingency_table(records)

    expect_true("N_ID" %in% names(result))
})

test_that("build_contingency_table N_ID sums to nrow of input", {
    records <- tibble::tibble(
        capture_1 = c(1L, 0L, 1L, 0L),
        capture_2 = c(0L, 1L, 0L, 1L)
    )
    result <- build_contingency_table(records)

    expect_equal(sum(result$N_ID), nrow(records))
})

test_that("build_contingency_table collapses duplicate rows", {
    records <- tibble::tibble(
        capture_1 = c(1L, 1L, 0L),
        capture_2 = c(0L, 0L, 1L)
    )
    result <- build_contingency_table(records)

    expect_lt(nrow(result), nrow(records))
    duplicated_row_count <- result$N_ID[
        result$capture_1 == 1L & result$capture_2 == 0L
    ]
    expect_equal(duplicated_row_count, 2L)
})

# ==============================================================================
# Group 5: validate_create_args() - input validation
# ==============================================================================

test_that("validate_create_args errors on non-numeric n_individuals", {
    expect_error(
        crcsim:::validate_create_args("10", 2, 1, FALSE),
        "n_individuals must be a numeric scalar"
    )
})

test_that("validate_create_args errors on non-scalar n_individuals", {
    expect_error(
        crcsim:::validate_create_args(c(10, 20), 2, 1, FALSE),
        "n_individuals must be a numeric scalar"
    )
})

test_that("validate_create_args errors on n_individuals less than 1", {
    expect_error(
        crcsim:::validate_create_args(0, 2, 1, FALSE),
        "n_individuals must be a positive integer"
    )
})

test_that("validate_create_args errors on non-integer n_individuals", {
    expect_error(
        crcsim:::validate_create_args(1.5, 2, 1, FALSE),
        "n_individuals must be a positive integer"
    )
})

test_that("validate_create_args errors on non-numeric n_captures", {
    expect_error(
        crcsim:::validate_create_args(10, "2", 1, FALSE),
        "n_captures must be a numeric scalar"
    )
})

test_that("validate_create_args errors on non-scalar n_captures", {
    expect_error(
        crcsim:::validate_create_args(10, c(2, 3), 1, FALSE),
        "n_captures must be a numeric scalar"
    )
})

test_that("validate_create_args errors on n_captures less than 1", {
    expect_error(
        crcsim:::validate_create_args(10, 0, 1, FALSE),
        "n_captures must be a positive integer"
    )
})

test_that("validate_create_args errors on non-integer n_captures", {
    expect_error(
        crcsim:::validate_create_args(10, 2.5, 1, FALSE),
        "n_captures must be a positive integer"
    )
})

# ==============================================================================
# Group 6: normalize_p_captures() - internal validation
# ==============================================================================

test_that("normalize_p_captures returns rep(0.5, n) when p_captures is NULL", {
    result <- crcsim:::normalize_p_captures(n_captures = 3, p_captures = NULL)

    expect_equal(result, rep(0.5, 3))
})

test_that("normalize_p_captures errors on non-numeric input", {
    expect_error(
        crcsim:::normalize_p_captures(2, c("0.3", "0.7")),
        "p_captures must be a numeric vector"
    )
})

test_that("normalize_p_captures errors when length does not match n_captures", {
    expect_error(
        crcsim:::normalize_p_captures(3, c(0.3, 0.7)),
        "Length of p_captures must equal n_captures"
    )
})

test_that("normalize_p_captures errors when any value is outside [0, 1]", {
    expect_error(
        crcsim:::normalize_p_captures(2, c(0.5, 1.5)),
        "All p_captures values must be between 0 and 1"
    )
})

test_that("normalize_p_captures returns valid input as numeric vector", {
    result <- crcsim:::normalize_p_captures(3, c(0.2, 0.5, 0.8))

    expect_equal(result, c(0.2, 0.5, 0.8))
    expect_type(result, "double")
})
