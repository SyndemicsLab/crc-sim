################################################################################
# File: test_loglinear.R                                                       #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

make_loglinear_fixture <- function(n_individuals = 300, n_captures = 3) {
    set.seed(123)
    sim_data <- create_data(
        n_individuals = n_individuals,
        n_captures = n_captures,
        return_frequency_table = TRUE
    )
    return(sim_data)
}

test_that("fit_loglinear_model returns a glm object with poisson", {
    data <- make_loglinear_fixture()
    formula_obj <- stats::as.formula("N_ID ~ capture_1 + capture_2")

    model <- fit_loglinear_model(
        data = data,
        formula_object = formula_obj,
        model_family = "poisson"
    )

    expect_s3_class(model, "glm")
})

test_that("fit_loglinear_model returns a negbin object with negbin", {
    data <- make_loglinear_fixture()
    formula_obj <- stats::as.formula("N_ID ~ capture_1 + capture_2")

    model <- fit_loglinear_model(
        data = data,
        formula_object = formula_obj,
        model_family = "negbin"
    )

    expect_s3_class(model, "glm")
})

test_that("fit_loglinear_model errors on invalid model_family", {
    data <- make_loglinear_fixture()
    formula_obj <- stats::as.formula("N_ID ~ capture_1 + capture_2")

    expect_error(
        fit_loglinear_model(
            data = data,
            formula_object = formula_obj,
            model_family = "invalid_family"
        ),
        "should be one of"
    )
})

test_that("fit_loglinear_model errors on invalid formula", {
    data <- make_loglinear_fixture()
    bad_formula <- stats::as.formula("N_ID ~ does_not_exist")

    expect_error(
        fit_loglinear_model(
            data = data,
            formula_object = bad_formula,
            model_family = "poisson"
        ),
        "does_not_exist"
    )
})
