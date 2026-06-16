################################################################################
# File: test_aic_selection.R                                                   #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-06-16                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

make_aic_fixture <- function(n_individuals = 300, n_captures = 3) {
    set.seed(123)
    sim_data <- create_data(
        n_individuals = n_individuals,
        n_captures = n_captures,
        return_frequency_table = TRUE
    )
    return(sim_data)
}

make_aic_options <- function(formulas = NULL, model = "poisson") {
    if (is.null(formulas)) {
        formulas <- list(
            stats::as.formula("N_ID ~ capture_1 + capture_2"),
            stats::as.formula("N_ID ~ capture_1 * capture_2")
        )
    }

    opts <- AICOptions$new(
        model = model,
        capture_columns = c("capture_1", "capture_2", "capture_3"),
        formula = formulas,
        frequency_col_name = "N_ID"
    )

    return(opts)
}

test_that("aic_selection errors when opts is not an AICOptions object", {
    model_data <- make_aic_fixture()

    expect_error(
        aic_selection(model_data, list()),
        "Invalid AICOptions object provided"
    )
})

test_that("aic_selection errors when data is not a frequency table", {
    set.seed(456)
    raw_data <- create_data(
        n_individuals = 100,
        n_captures = 3,
        return_frequency_table = FALSE
    )

    opts <- make_aic_options()

    expect_error(
        aic_selection(raw_data, opts),
        "Data must be a frequency table"
    )
})

test_that("aic_selection errors when frequency column exists but is not numeric", {
    model_data <- make_aic_fixture()
    model_data$N_ID <- as.character(model_data$N_ID)
    opts <- make_aic_options()

    expect_error(
        aic_selection(model_data, opts),
        "Data must be a frequency table"
    )
})

test_that("build_aic_error_row returns expected NA fields and error message", {
    formula_object <- stats::as.formula("N_ID ~ capture_1 + capture_2")
    err <- simpleError("synthetic model fitting error")

    out <- build_aic_error_row(formula_object, err)

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 1)
    expect_named(
        out,
        c(
            "formula",
            "estimate",
            "AIC",
            "lower_ci",
            "upper_ci",
            "error"
        )
    )
    expect_true(is.na(out$estimate))
    expect_true(is.na(out$AIC))
    expect_true(is.na(out$lower_ci))
    expect_true(is.na(out$upper_ci))
    expect_equal(out$error, "synthetic model fitting error")
})

test_that("build_aic_success_row returns expected schema and rounded values", {
    model_data <- make_aic_fixture()
    formula_object <- stats::as.formula("N_ID ~ capture_1 + capture_2")

    model <- fit_loglinear_model(
        model_data = model_data,
        formula_object = formula_object,
        model_family = "poisson"
    )

    out <- build_aic_success_row(model, formula_object)

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 1)
    expect_named(
        out,
        c(
            "formula",
            "estimate",
            "AIC",
            "lower_ci",
            "upper_ci",
            "error"
        )
    )
    expect_type(out$estimate, "double")
    expect_type(out$AIC, "double")
    expect_true(is.na(out$error))
    expect_lte(out$lower_ci, out$estimate)
    expect_lte(out$estimate, out$upper_ci)
})

test_that("evaluate_formula_with_aic returns error-row shape on failures", {
    model_data <- make_aic_fixture()
    bad_formula <- stats::as.formula("N_ID ~ does_not_exist")

    out <- evaluate_formula_with_aic(
        formula_object = bad_formula,
        model_data = model_data,
        model_family = "poisson"
    )

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 1)
    expect_true(is.na(out$estimate))
    expect_true(is.na(out$AIC))
    expect_match(out$error, "does_not_exist")
})

test_that("evaluate_formula_with_aic rejects unsupported model family", {
    model_data <- make_aic_fixture()
    formula_object <- stats::as.formula("N_ID ~ capture_1")

    expect_error(
        evaluate_formula_with_aic(
            formula_object = formula_object,
            model_data = model_data,
            model_family = "badfamily"
        ),
        "should be one of"
    )
})

test_that(
    paste(
        "aic_selection returns one row per formula with required columns",
        "and expected field types"
    ),
    {
        model_data <- make_aic_fixture()
        formulas <- list(
            stats::as.formula("N_ID ~ capture_1 + capture_2"),
            stats::as.formula("N_ID ~ capture_1 * capture_2")
        )
        opts <- make_aic_options(formulas = formulas)

        out <- aic_selection(model_data, opts)

        expect_s3_class(out, "data.frame")
        expect_equal(nrow(out), length(formulas))
        expect_named(
            out,
            c(
                "formula",
                "estimate",
                "AIC",
                "lower_ci",
                "upper_ci",
                "error"
            )
        )
        expect_true(all(vapply(out$formula, is.character, logical(1))))
    }
)

test_that("aic_selection handles mixed success and error formulas", {
    model_data <- make_aic_fixture()
    formulas <- list(
        stats::as.formula("N_ID ~ capture_1 + capture_2"),
        stats::as.formula("N_ID ~ does_not_exist")
    )
    opts <- make_aic_options(formulas = formulas)

    out <- aic_selection(model_data, opts)

    expect_equal(nrow(out), 2)
    expect_true(anyNA(out$error))
    expect_false(all(is.na(out$error)))
})

test_that("aic_selection output is sorted by ascending AIC for successful rows", {
    model_data <- make_aic_fixture()
    formulas <- list(
        stats::as.formula("N_ID ~ capture_1 + capture_2"),
        stats::as.formula("N_ID ~ capture_1 + capture_2 + capture_3"),
        stats::as.formula("N_ID ~ capture_1 * capture_2")
    )
    opts <- make_aic_options(formulas = formulas)

    out <- aic_selection(model_data, opts)

    ok <- out[!is.na(out$AIC), , drop = FALSE]
    if (nrow(ok) > 1) {
        expect_true(all(diff(ok$AIC) >= 0))
    } else {
        succeed()
    }
})

test_that("aic_selection result contains formula strings from provided formulas", {
    model_data <- make_aic_fixture()
    formula_1 <- stats::as.formula("N_ID ~ capture_1 + capture_2")
    formula_2 <- stats::as.formula("N_ID ~ capture_1 * capture_2")
    opts <- make_aic_options(formulas = list(formula_1, formula_2))

    out <- aic_selection(model_data, opts)

    expect_true(any(grepl("capture_1", out$formula, fixed = TRUE)))
    expect_true(any(grepl("capture_2", out$formula, fixed = TRUE)))
})

test_that("aic_selection preserves error text in output rows", {
    model_data <- make_aic_fixture()
    formulas <- list(stats::as.formula("N_ID ~ not_a_real_column"))
    opts <- make_aic_options(formulas = formulas)

    out <- aic_selection(model_data, opts)

    expect_equal(nrow(out), 1)
    expect_true(is.na(out$estimate))
    expect_true(is.na(out$AIC))
    expect_false(is.na(out$error))
    expect_match(out$error, "not_a_real_column")
})
