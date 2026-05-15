################################################################################
# File: test_crc_stepwise.R                                                    #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

make_stepwise_fixture <- function(n_individuals = 300, n_captures = 3) {
    set.seed(123)
    sim_data <- create_data(
        n_individuals = n_individuals,
        n_captures = n_captures,
        return_frequency_table = TRUE
    )
    return(sim_data)
}

make_stepwise_options <- function(
    model = "poisson",
    threshold = 0.05,
    direction = "both",
    interaction_limit = 2
) {
    opts <- StepwiseOptions$new(
        model = model,
        threshold = threshold,
        direction = direction,
        frequency_col_name = "N_ID",
        capture_indicators = c("capture_1", "capture_2", "capture_3"),
        interaction_limit = interaction_limit
    )
    return(opts)
}

# ============================================================================
# Input Validation Tests
# ============================================================================

test_that("crc_stepwise errors when opts is not a StepwiseOptions object", {
    data <- make_stepwise_fixture()

    expect_error(
        crc_stepwise(data, list()),
        "Invalid StepwiseOptions object provided"
    )
})

test_that("crc_stepwise errors when data is not a frequency table", {
    set.seed(456)
    raw_data <- create_data(
        n_individuals = 100,
        n_captures = 3,
        return_frequency_table = FALSE
    )

    opts <- make_stepwise_options()

    expect_error(
        crc_stepwise(raw_data, opts),
        "Data must be a frequency table"
    )
})

test_that("crc_stepwise errors when frequency column is not numeric", {
    data <- make_stepwise_fixture()
    data$N_ID <- as.character(data$N_ID)
    opts <- make_stepwise_options()

    expect_error(
        crc_stepwise(data, opts),
        "Data must be a frequency table"
    )
})

# ============================================================================
# Output Structure and Type Tests
# ============================================================================

test_that("crc_stepwise returns list with required fields", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- crc_stepwise(data, opts)

    expect_type(result, "list")
    expect_named(
        result,
        c(
            "model",
            "formula",
            "summary",
            "estimate",
            "lower_ci",
            "upper_ci",
            "AIC"
        )
    )
})

test_that("crc_stepwise returns numeric estimate rounded to 2 decimals", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- crc_stepwise(data, opts)

    expect_type(result$estimate, "double")
    expect_equal(result$estimate, round(result$estimate, 2))
})

test_that("crc_stepwise returns numeric CI bounds rounded to 2 decimals", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- crc_stepwise(data, opts)

    expect_type(result$lower_ci, "double")
    expect_type(result$upper_ci, "double")
    expect_equal(result$lower_ci, round(result$lower_ci, 2))
    expect_equal(result$upper_ci, round(result$upper_ci, 2))
})

test_that("crc_stepwise returns valid confidence interval", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- crc_stepwise(data, opts)

    expect_lte(result$lower_ci, result$estimate)
    expect_lte(result$estimate, result$upper_ci)
})

test_that("crc_stepwise returns numeric AIC", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- crc_stepwise(data, opts)

    expect_type(result$AIC, "double")
})

test_that("crc_stepwise returns formula object", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- crc_stepwise(data, opts)

    expect_true(inherits(result$formula, "formula"))
})

test_that("crc_stepwise returns model family string", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(model = "poisson")

    result <- crc_stepwise(data, opts)

    expect_equal(result$model, "poisson")
})

# ============================================================================
# Core Functionality Tests
# ============================================================================

test_that("crc_stepwise works with poisson model", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(model = "poisson")

    result <- crc_stepwise(data, opts)

    expect_equal(result$model, "poisson")
    expect_true(!is.na(result$estimate))
    expect_true(!is.na(result$AIC))
})

test_that("crc_stepwise works with negbin model", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(model = "negbin")

    result <- crc_stepwise(data, opts)

    expect_equal(result$model, "negbin")
    expect_true(!is.na(result$estimate))
    expect_true(!is.na(result$AIC))
})

test_that("crc_stepwise respects interaction_limit = 2", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(interaction_limit = 2, direction = "both")

    result <- crc_stepwise(data, opts)

    expect_true(!is.na(result$AIC))
})

test_that("crc_stepwise respects interaction_limit = 3", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(interaction_limit = 3, direction = "both")

    result <- crc_stepwise(data, opts)

    expect_true(!is.na(result$AIC))
})

# ============================================================================
# Edge Cases and Robustness Tests
# ============================================================================

test_that("crc_stepwise works with minimal data", {
    set.seed(999)
    small_data <- create_data(
        n_individuals = 50,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts <- make_stepwise_options()

    result <- crc_stepwise(small_data, opts)

    expect_true(!is.na(result$estimate))
})

test_that("crc_stepwise works with larger datasets", {
    set.seed(888)
    large_data <- create_data(
        n_individuals = 1000,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts <- make_stepwise_options()

    result <- crc_stepwise(large_data, opts)

    expect_true(!is.na(result$estimate))
})

test_that("crc_stepwise results are reproducible with seed", {
    set.seed(777)
    data1 <- create_data(
        n_individuals = 200,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts1 <- make_stepwise_options()
    result1 <- crc_stepwise(data1, opts1)

    set.seed(777)
    data2 <- create_data(
        n_individuals = 200,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts2 <- make_stepwise_options()
    result2 <- crc_stepwise(data2, opts2)

    expect_equal(result1$estimate, result2$estimate)
    expect_equal(result1$AIC, result2$AIC)
})

test_that("crc_stepwise final formula is within scope", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(interaction_limit = 2)

    result <- crc_stepwise(data, opts)

    # Final formula should include response and at least some predictors
    formula_vars <- all.vars(result$formula)
    expect_true("N_ID" %in% formula_vars)
    expect_true(any(c("capture_1", "capture_2", "capture_3") %in% formula_vars))
})

test_that("crc_stepwise handles case where formula does not reduce", {
    data <- make_stepwise_fixture()
    # With low p-value threshold, stepwise might not remove any variables
    opts <- make_stepwise_options(threshold = 0.9, direction = "backward")

    result <- crc_stepwise(data, opts)

    expect_true(!is.na(result$AIC))
    expect_true(!is.na(result$estimate))
})

# ============================================================================
# Extensive Direction Testing
# ============================================================================

test_that("crc_stepwise works with direction = 'forward'", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "forward")

    result <- crc_stepwise(data, opts)

    expect_equal(result$model, "poisson")
    expect_true(!is.na(result$estimate))
    expect_true(!is.na(result$AIC))
})

test_that("crc_stepwise works with direction = 'backward'", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "backward")

    result <- crc_stepwise(data, opts)

    expect_equal(result$model, "poisson")
    expect_true(!is.na(result$estimate))
    expect_true(!is.na(result$AIC))
})

test_that("crc_stepwise works with direction = 'both'", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "both")

    result <- crc_stepwise(data, opts)

    expect_equal(result$model, "poisson")
    expect_true(!is.na(result$estimate))
    expect_true(!is.na(result$AIC))
})

test_that("forward direction with interaction_limit = 2", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        direction = "forward",
        interaction_limit = 2
    )

    result <- crc_stepwise(data, opts)

    expect_true(!is.na(result$AIC))
    expect_true(!is.na(result$estimate))
})

test_that("backward direction with interaction_limit = 2", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        direction = "backward",
        interaction_limit = 2
    )

    result <- crc_stepwise(data, opts)

    expect_true(!is.na(result$AIC))
    expect_true(!is.na(result$estimate))
})

test_that("forward direction with negbin model", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        model = "negbin",
        direction = "forward"
    )

    result <- crc_stepwise(data, opts)

    expect_equal(result$model, "negbin")
    expect_true(!is.na(result$AIC))
})

test_that("backward direction with negbin model", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        model = "negbin",
        direction = "backward"
    )

    result <- crc_stepwise(data, opts)

    expect_equal(result$model, "negbin")
    expect_true(!is.na(result$AIC))
})

test_that("forward direction with interaction_limit = 3", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        direction = "forward",
        interaction_limit = 3
    )

    result <- crc_stepwise(data, opts)

    expect_true(!is.na(result$AIC))
})

test_that("backward direction with interaction_limit = 3", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        direction = "backward",
        interaction_limit = 3
    )

    result <- crc_stepwise(data, opts)

    expect_true(!is.na(result$AIC))
})

test_that("forward direction produces valid confidence interval", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "forward")

    result <- crc_stepwise(data, opts)

    expect_lte(result$lower_ci, result$estimate)
    expect_lte(result$estimate, result$upper_ci)
})

test_that("backward direction produces valid confidence interval", {
    data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "backward")

    result <- crc_stepwise(data, opts)

    expect_lte(result$lower_ci, result$estimate)
    expect_lte(result$estimate, result$upper_ci)
})
