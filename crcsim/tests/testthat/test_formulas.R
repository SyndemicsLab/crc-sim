################################################################################
# File: test_formulas.R                                                        #
# Project: crcsim                                                              #
# Created Date: 2026-08-31                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("is_formula identifies formula objects", {
    expect_true(crcsim:::is_formula(stats::as.formula("y ~ x")))
    expect_false(crcsim:::is_formula("y ~ x"))
    expect_false(crcsim:::is_formula(NULL))
    expect_false(crcsim:::is_formula(1))
})

test_that("formula_list generates all two-predictor formula variants", {
    result <- crcsim:::formula_list("y", c("a", "b"))

    expect_type(result, "list")
    expect_length(result, 6)
    expect_true(all(vapply(result, crcsim:::is_formula, logical(1))))
    expect_identical(
        vapply(result, deparse, character(1)),
        c(
            "y ~ a",
            "y ~ b",
            "y ~ a + b",
            "y ~ a * b",
            "y ~ a + a * b",
            "y ~ b + a * b"
        )
    )
})

test_that("formula_list generates unique higher-order formulas", {
    result <- crcsim:::formula_list("outcome", c("x1", "x2", "x3"))
    formula_text <- vapply(result, deparse, character(1))

    expect_length(result, 23)
    expect_length(unique(formula_text), length(formula_text))
    expect_true(any(grepl("x1 \\* x2 \\* x3", formula_text)))
    expect_true(any(grepl("outcome ~ x1 \\+ x2 \\+ x3", formula_text)))
})

test_that("formula_list errors for fewer than two predictors", {
    expect_error(crcsim:::formula_list("y", character()), "n < m")
    expect_error(crcsim:::formula_list("y", "x"), "n < m")
})
