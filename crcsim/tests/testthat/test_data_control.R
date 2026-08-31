################################################################################
# File: test_data_control_unit.R                                               #
# Project: crcsim                                                              #
# Created Date: 2026-08-28                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("sample_single_capture returns one binary value per probability", {
    set.seed(1)
    result <- crcsim:::sample_single_capture(c(0.2, 0.5, 0.8))

    expect_type(result, "integer")
    expect_length(result, 3)
    expect_true(all(result %in% c(0L, 1L)))
})

test_that("sample_single_capture respects deterministic probabilities", {
    result <- crcsim:::sample_single_capture(c(0, 1, 0, 1))

    expect_identical(result, c(0L, 1L, 0L, 1L))
})

test_that("sample_captures returns named binary capture columns", {
    result <- sample_captures(4, c(0, 1, 0.5))

    expect_s3_class(result, "data.frame")
    expect_equal(dim(result), c(4, 3))
    expect_named(result, c("capture_1", "capture_2", "capture_3"))
    expect_true(all(result$capture_1 == 0L))
    expect_true(all(result$capture_2 == 1L))
    expect_true(all(unlist(result$capture_3) %in% c(0L, 1L)))
})

test_that("sample_captures validates capture probabilities through create validation", {
    expect_error(sample_captures(3, numeric()), "positive integer")
    expect_error(sample_captures(3, c(0.2, 1.2)), "between 0 and 1")
})

test_that("captures_from_covariates adds requested capture columns", {
    data <- data.frame(covariate_1 = 1:4)
    result <- captures_from_covariates(data, 2, function(x) rep(1, nrow(x)))

    expect_named(result, c("covariate_1", "capture_1", "capture_2"))
    expect_true(all(result$capture_1 == 1L))
    expect_true(all(result$capture_2 == 1L))
})

test_that("captures_from_covariates does not overwrite existing captures", {
    data <- data.frame(capture_1 = c(0L, 1L), covariate_1 = c(2, 3))
    result <- captures_from_covariates(
        data,
        2,
        function(x) rep(1, nrow(x))
    )

    expect_identical(result$capture_1, data$capture_1)
    expect_true(all(result$capture_2 == 1L))
})

test_that("captures_from_covariates returns input for no requested captures", {
    data <- data.frame(covariate_1 = 1:2)

    expect_identical(
        captures_from_covariates(data, 0, function(x) stop("not called")),
        data
    )
})

test_that("captures_from_covariates validates inputs and probabilities", {
    data <- data.frame(covariate_1 = 1:2)

    expect_error(
        captures_from_covariates(1:2, 1, function(x) 0.5),
        "data_table must be a data.frame"
    )
    expect_error(
        captures_from_covariates(data, 1.5, function(x) 0.5),
        "n_cols must be an integer scalar"
    )
    expect_error(
        captures_from_covariates(data, 1, 0.5),
        "cov_func must be a function"
    )
    expect_error(
        captures_from_covariates(data, 1, function(x) 0.5),
        "length nrow"
    )
    expect_error(
        captures_from_covariates(data, 1, function(x) c(0.5, NA)),
        "between 0 and 1"
    )
})

test_that("covariate_conditions returns probabilities in the unit interval", {
    data <- data.frame(
        covariate_1 = c(0, 1),
        covariate_2 = c(1, 0),
        covariate_3 = c(0, 1),
        covariate_4 = c(1, 0)
    )

    result <- covariate_conditions(data)

    expect_type(result, "double")
    expect_length(result, 2)
    expect_true(all(result > 0 & result < 1))
})

test_that("suppress_data removes frequencies at or below threshold", {
    data <- data.frame(id = 1:4, count = c(1, 2, 3, 4))

    result <- suppress_data(data, 2, freq_col = "count")

    expect_equal(result$id, c(3L, 4L))
    expect_equal(result$count, c(3, 4))
})

test_that("suppress_data validates its inputs", {
    data <- data.frame(N_ID = 1:2)

    expect_error(suppress_data(1:2, 1), "data must be a data.frame")
    expect_error(suppress_data(data, c(1, 2)), "numeric scalar")
    expect_error(suppress_data(data, 1, c("N_ID", "other")), "single column")
    expect_error(
        suppress_data(data, 1, "missing"),
        "Frequency column not found"
    )
    expect_error(
        suppress_data(data.frame(N_ID = c("1", "2")), 1),
        "Frequency column must be numeric"
    )
})

test_that("capture row extraction separates captured and uncaptured rows", {
    data <- tibble::tibble(
        capture_1 = c(0, 1, 0, 1),
        capture_2 = c(0, 0, 1, 1),
        value = 1:4
    )

    captured <- expect_no_warning(
        extract_captured_data(data, c("capture_1", "capture_2"))
    )
    uncaptured <- expect_no_warning(
        extract_uncaptured_data(data, c("capture_1", "capture_2"))
    )

    expect_equal(captured$value, c(2, 3, 4))
    expect_equal(uncaptured$value, 1)
    expect_s3_class(captured, "tbl_df")
    expect_s3_class(uncaptured, "tbl_df")
})

test_that("build_contingency_table counts duplicate records", {
    records <- tibble::tibble(
        capture_1 = c(0L, 0L, 1L),
        capture_2 = c(1L, 1L, 0L),
        covariate_1 = c(2, 2, 3)
    )

    result <- build_contingency_table(records)

    expect_named(
        result,
        c("capture_1", "capture_2", "covariate_1", "N_ID")
    )
    expect_equal(sum(result$N_ID), nrow(records))
    expect_equal(nrow(result), 2)
    expect_equal(
        result$N_ID[result$capture_1 == 0L & result$capture_2 == 1L],
        2L
    )
})

test_that("probability predicates accept valid values and reject invalid ones", {
    expect_true(is_probability_vector(c(0, 0.5, 1)))
    expect_false(is_probability_vector(c(0, 1.1)))
    expect_false(is_probability_vector(c("0", "1")))
    expect_true(is_probability_scalar(0.5))
    expect_false(is_probability_scalar(c(0.2, 0.8)))
    expect_false(is_probability_scalar(-0.1))
})

test_that("data predicates identify binary columns and frequency tables", {
    data <- data.frame(binary = c(0, 1), count = c(2, 3))

    expect_true(is_column_binary(data$binary))
    expect_false(is_column_binary(c(0, 0)))
    expect_true(is_frequency_table(data, "count"))
    expect_false(is_frequency_table(data, "missing"))
    expect_false(is_frequency_table(
        transform(data, count = as.character(count)),
        "count"
    ))
})

test_that("all_int_cols_to_numeric converts only integer columns", {
    data <- data.frame(integer_col = 1:2, numeric_col = c(1.5, 2.5))

    result <- all_int_cols_to_numeric(data)

    expect_type(result$integer_col, "double")
    expect_type(result$numeric_col, "double")
    expect_equal(result$integer_col, c(1, 2))
    expect_equal(result$numeric_col, c(1.5, 2.5))
})

test_that("create_data generates deterministic boundary captures", {
    result <- create_data(
        n_individuals = 5,
        n_captures = 2,
        p_captures = c(0, 1)
    )

    expect_s3_class(result, "tbl_df")
    expect_equal(nrow(result), 5)
    expect_true(all(result$capture_1 == 0L))
    expect_true(all(result$capture_2 == 1L))
})

test_that("create_data generates configured covariate columns and types", {
    specs <- data.frame(
        distribution = c("uniform", "normal"),
        p1 = c(0, 10),
        p2 = c(1, 1),
        dtype = c("integer", "numeric")
    )

    set.seed(1)
    result <- create_data(20, 1, covariate_ranges = specs)

    expect_named(
        result,
        c("capture_1", "covariate_1", "covariate_2")
    )
    expect_type(result$covariate_1, "integer")
    expect_type(result$covariate_2, "double")
    expect_true(all(result$covariate_1 >= 0 & result$covariate_1 <= 1))
})

test_that("create_data frequency output preserves total population", {
    result <- create_data(100, 2, return_frequency_table = TRUE)

    expect_s3_class(result, "tbl_df")
    expect_true("N_ID" %in% names(result))
    expect_equal(sum(result$N_ID), 100)
})

test_that("create_data validates top-level arguments", {
    expect_error(create_data("10", 2), "n_individuals must be a numeric scalar")
    expect_error(create_data(1.5, 2), "positive integer")
    expect_error(create_data(10, 0), "positive integer")
    expect_error(
        create_data(10, 2, covariate_ranges = 1:2),
        "data.frame or NULL"
    )
    expect_error(
        create_data(10, 2, return_frequency_table = 1),
        "logical scalar"
    )
})

test_that("covariate specification normalization validates and standardizes input", {
    specs <- data.frame(
        distribution = c("UNIFORM", "NORMAL"),
        p1 = c(0, 2),
        p2 = c(1, 0.5),
        dtype = c("INTEGER", "NUMERIC"),
        ignored = 1:2
    )

    result <- crcsim:::normalize_covariate_ranges(specs)

    expect_named(result, c("distribution", "p1", "p2", "dtype"))
    expect_equal(result$distribution, c("uniform", "normal"))
    expect_equal(result$dtype, c("integer", "numeric"))
})

test_that("covariate specification normalization handles null and invalid input", {
    empty <- crcsim:::normalize_covariate_ranges(NULL)

    expect_equal(names(empty), c("distribution", "p1", "p2", "dtype"))
    expect_equal(nrow(empty), 0)
    expect_error(
        crcsim:::normalize_covariate_ranges(data.frame(p1 = 0, p2 = 1)),
        "must contain columns"
    )
    expect_error(
        crcsim:::normalize_covariate_ranges(data.frame(
            distribution = "poisson",
            p1 = 0,
            p2 = 1,
            dtype = "numeric"
        )),
        "distribution must be one of"
    )
    expect_error(
        crcsim:::normalize_covariate_ranges(data.frame(
            distribution = "uniform",
            p1 = 2,
            p2 = 1,
            dtype = "numeric"
        )),
        "p1 \\(low\\) must be <= p2"
    )
    expect_error(
        crcsim:::normalize_covariate_ranges(data.frame(
            distribution = "normal",
            p1 = 0,
            p2 = 0,
            dtype = "numeric"
        )),
        "standard deviation\\) must be > 0"
    )
})

test_that("sample_covariate_value supports uniform and normal output types", {
    uniform <- data.frame(
        distribution = "uniform",
        p1 = 2,
        p2 = 3,
        dtype = "integer"
    )
    normal <- data.frame(
        distribution = "normal",
        p1 = 10,
        p2 = 0.1,
        dtype = "numeric"
    )

    set.seed(1)
    uniform_value <- crcsim:::sample_covariate_value(uniform)
    normal_value <- crcsim:::sample_covariate_value(normal)

    expect_type(uniform_value, "integer")
    expect_true(uniform_value %in% 2:3)
    expect_type(normal_value, "double")
    expect_true(is.finite(normal_value))
})
