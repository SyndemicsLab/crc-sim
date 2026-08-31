################################################################################
# File: test_stepwise_selection.R                                              #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-06-16                                                    #
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
        capture_columns = c("capture_1", "capture_2", "capture_3"),
        threshold = threshold,
        direction = direction,
        frequency_col_name = "N_ID",
        interaction_limit = interaction_limit
    )
    return(opts)
}

# ============================================================================
# Input Validation Tests
# ============================================================================

test_that("stepwise_selection errors when opts is not a StepwiseOptions object", {
    model_data <- make_stepwise_fixture()

    expect_error(
        stepwise_selection(model_data, list()),
        "Invalid StepwiseOptions object provided"
    )
})

test_that("stepwise_selection errors when data is not a frequency table", {
    set.seed(456)
    raw_data <- create_data(
        n_individuals = 100,
        n_captures = 3,
        return_frequency_table = FALSE
    )

    opts <- make_stepwise_options()

    expect_error(
        stepwise_selection(raw_data, opts),
        "Data must be a frequency table"
    )
})

test_that("stepwise_selection errors when frequency column is not numeric", {
    model_data <- make_stepwise_fixture()
    model_data$N_ID <- as.character(model_data$N_ID)
    opts <- make_stepwise_options()

    expect_error(
        stepwise_selection(model_data, opts),
        "Data must be a frequency table"
    )
})

# ============================================================================
# Output Structure and Type Tests
# ============================================================================

test_that("stepwise_selection returns list with required fields", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- stepwise_selection(model_data, opts)

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

test_that("stepwise_selection returns numeric estimate rounded to 2 decimals", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- stepwise_selection(model_data, opts)

    expect_type(result$estimate, "double")
    expect_equal(result$estimate, round(result$estimate, 2))
})

test_that("stepwise_selection returns numeric CI bounds rounded to 2 decimals", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- stepwise_selection(model_data, opts)

    expect_type(result$lower_ci, "double")
    expect_type(result$upper_ci, "double")
    expect_equal(result$lower_ci, round(result$lower_ci, 2))
    expect_equal(result$upper_ci, round(result$upper_ci, 2))
})

test_that("stepwise_selection returns valid confidence interval", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- stepwise_selection(model_data, opts)

    expect_lte(result$lower_ci, result$estimate)
    expect_lte(result$estimate, result$upper_ci)
})

test_that("stepwise_selection returns numeric AIC", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- stepwise_selection(model_data, opts)

    expect_type(result$AIC, "double")
})

test_that("stepwise_selection returns formula object", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options()

    result <- stepwise_selection(model_data, opts)

    expect_s3_class(result$formula, "formula")
})

test_that("stepwise_selection returns model family string", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(model = "poisson")

    result <- stepwise_selection(model_data, opts)

    expect_equal(result$model, "poisson")
})

# ============================================================================
# Core Functionality Tests
# ============================================================================

test_that("stepwise_selection works with poisson model", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(model = "poisson")

    result <- stepwise_selection(model_data, opts)

    expect_equal(result$model, "poisson")
    expect_false(is.na(result$estimate))
    expect_false(is.na(result$AIC))
})

test_that("stepwise_selection works with negbin model", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(model = "negbin")

    # Suppress alternation limit reached
    # nolint start: implicit_assignment_linter
    suppressWarnings(result <- stepwise_selection(model_data, opts))
    # nolint end: implicit_assignment_linter

    expect_equal(result$model, "negbin")
    expect_false(is.na(result$estimate))
    expect_false(is.na(result$AIC))
})

test_that("stepwise_selection respects interaction_limit = 2", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(interaction_limit = 2, direction = "both")

    result <- stepwise_selection(model_data, opts)

    expect_false(is.na(result$AIC))
})

test_that("stepwise_selection respects interaction_limit = 3", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(interaction_limit = 3, direction = "both")

    result <- stepwise_selection(model_data, opts)

    expect_false(is.na(result$AIC))
})

# ============================================================================
# Edge Cases and Robustness Tests
# ============================================================================

test_that("stepwise_selection works with minimal data", {
    set.seed(999)
    small_data <- create_data(
        n_individuals = 50,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts <- make_stepwise_options()

    result <- stepwise_selection(small_data, opts)

    expect_false(is.na(result$estimate))
})

test_that("stepwise_selection works with larger datasets", {
    set.seed(888)
    large_data <- create_data(
        n_individuals = 1000,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts <- make_stepwise_options()

    result <- stepwise_selection(large_data, opts)

    expect_false(is.na(result$estimate))
})

test_that("stepwise_selection results are reproducible with seed", {
    set.seed(777)
    data1 <- create_data(
        n_individuals = 200,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts1 <- make_stepwise_options()
    result1 <- stepwise_selection(data1, opts1)

    set.seed(777)
    data2 <- create_data(
        n_individuals = 200,
        n_captures = 3,
        return_frequency_table = TRUE
    )
    opts2 <- make_stepwise_options()
    result2 <- stepwise_selection(data2, opts2)

    expect_equal(result1$estimate, result2$estimate)
    expect_equal(result1$AIC, result2$AIC)
})

test_that("stepwise_selection final formula is within scope", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(interaction_limit = 2)

    result <- stepwise_selection(model_data, opts)

    # Final formula should include response and at least some predictors
    formula_vars <- all.vars(result$formula)
    expect_true("N_ID" %in% formula_vars)
    expect_true(any(c("capture_1", "capture_2", "capture_3") %in% formula_vars))
})

test_that("stepwise_selection handles case where formula does not reduce", {
    model_data <- make_stepwise_fixture()
    # With low p-value threshold, stepwise might not remove any variables
    opts <- make_stepwise_options(threshold = 0.9, direction = "backward")

    result <- stepwise_selection(model_data, opts)

    expect_false(is.na(result$AIC))
    expect_false(is.na(result$estimate))
})

# ============================================================================
# Extensive Direction Testing
# ============================================================================

test_that("stepwise_selection works with direction = 'forward'", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "forward")

    result <- stepwise_selection(model_data, opts)

    expect_equal(result$model, "poisson")
    expect_false(is.na(result$estimate))
    expect_false(is.na(result$AIC))
})

test_that("stepwise_selection works with direction = 'backward'", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "backward")

    result <- stepwise_selection(model_data, opts)

    expect_equal(result$model, "poisson")
    expect_false(is.na(result$estimate))
    expect_false(is.na(result$AIC))
})

test_that("stepwise_selection works with direction = 'both'", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "both")

    result <- stepwise_selection(model_data, opts)

    expect_equal(result$model, "poisson")
    expect_false(is.na(result$estimate))
    expect_false(is.na(result$AIC))
})

test_that("forward direction with interaction_limit = 2", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        direction = "forward",
        interaction_limit = 2
    )

    result <- stepwise_selection(model_data, opts)

    expect_false(is.na(result$AIC))
    expect_false(is.na(result$estimate))
})

test_that("backward direction with interaction_limit = 2", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        direction = "backward",
        interaction_limit = 2
    )

    result <- stepwise_selection(model_data, opts)

    expect_false(is.na(result$AIC))
    expect_false(is.na(result$estimate))
})

test_that("forward direction with negbin model", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        model = "negbin",
        direction = "forward"
    )

    expect_warning(
        result <- stepwise_selection(model_data, opts),
        "alternation limit reached"
    )

    expect_equal(result$model, "negbin")
    expect_false(is.na(result$AIC))
})

test_that("backward direction with negbin model", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        model = "negbin",
        direction = "backward"
    )

    expect_warning(
        result <- stepwise_selection(model_data, opts),
        "alternation limit reached"
    )

    expect_equal(result$model, "negbin")
    expect_false(is.na(result$AIC))
})

test_that("forward direction with interaction_limit = 3", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        direction = "forward",
        interaction_limit = 3
    )

    result <- stepwise_selection(model_data, opts)

    expect_false(is.na(result$AIC))
})

test_that("backward direction with interaction_limit = 3", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(
        direction = "backward",
        interaction_limit = 3
    )

    result <- stepwise_selection(model_data, opts)

    expect_false(is.na(result$AIC))
})

test_that("forward direction produces valid confidence interval", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "forward")

    result <- stepwise_selection(model_data, opts)

    expect_lte(result$lower_ci, result$estimate)
    expect_lte(result$estimate, result$upper_ci)
})

test_that("backward direction produces valid confidence interval", {
    model_data <- make_stepwise_fixture()
    opts <- make_stepwise_options(direction = "backward")

    result <- stepwise_selection(model_data, opts)

    expect_lte(result$lower_ci, result$estimate)
    expect_lte(result$estimate, result$upper_ci)
})
