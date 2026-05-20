################################################################################
# File: test_options.R                                                         #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-19                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("AICOptions returns the provided formula when formula is present", {
    supplied_formula <- stats::as.formula("N_ID ~ capture_1 + capture_2")

    opts <- AICOptions$new(
        model = "poisson",
        frequency_col_name = "N_ID",
        capture_indicators = c("capture_1", "capture_2"),
        formula = supplied_formula
    )

    expect_equal(opts$formulas, supplied_formula)
})

test_that("AICOptions creates formulas when formula is absent", {
    expected_formulas <- formula_list("N_ID", c("capture_1", "capture_2"))

    opts <- AICOptions$new(
        model = "poisson",
        frequency_col_name = "N_ID",
        capture_indicators = c("capture_1", "capture_2")
    )

    expect_equal(opts$formulas, expected_formulas)
})

test_that(
    paste(
        "AICOptions errors when formula is absent and required inputs",
        "are missing"
    ),
    {
        expect_error(
            AICOptions$new(model = "poisson"),
            paste(
                "If formula is not provided, frequency_col_name and",
                "binary_variables must be specified."
            )
        )
    }
)

test_that("AICOptions prioritizes formula over generated formulas", {
    supplied_formula <- stats::as.formula("N_ID ~ capture_1 * capture_2")
    generated_formulas <- formula_list("N_ID", c("capture_1", "capture_2"))

    opts <- AICOptions$new(
        model = "poisson",
        frequency_col_name = "N_ID",
        capture_indicators = c("capture_1", "capture_2"),
        formula = supplied_formula
    )

    expect_equal(opts$formulas, supplied_formula)
    expect_false(identical(opts$formulas, generated_formulas))
})

# Test 1: Options constructor smoke test
test_that(
    paste(
        "Options constructor stores model and threshold fields",
        "exactly as provided"
    ),
    {
        opts <- Options$new(model = "poisson", threshold = 0.05)

        expect_equal(opts$model, "poisson")
        expect_equal(opts$threshold, 0.05)
    }
)

# Test 2: FrequencyOptions constructor and mutators
test_that(
    paste(
        "FrequencyOptions initializes with NULL formulas",
        "and the add_formula and add_formulas mutators evolve formulas"
    ),
    {
        opts <- FrequencyOptions$new(
            model = "poisson",
            threshold = 0.05,
            formulas = NULL,
            frequency_col_name = "N_ID"
        )

        expect_null(opts$formulas)

        # Add first formula
        formula1 <- stats::as.formula("N_ID ~ capture_1")
        result1 <- opts$add_formula(formula1)
        expect_equal(opts$formulas, formula1)

        # Add second formula
        formula2 <- stats::as.formula("N_ID ~ capture_2")
        opts$add_formula(formula2)
        expect_length(opts$formulas, 2)

        # Add multiple formulas
        formula3 <- stats::as.formula("N_ID ~ capture_3")
        formula4 <- stats::as.formula("N_ID ~ capture_4")
        opts$add_formulas(c(formula3, formula4))
        expect_length(opts$formulas, 4)
    }
)

# Test 3: StepwiseOptions constructor mapping
test_that(
    paste(
        "StepwiseOptions initializes all fields and preserves",
        "inherited model and threshold fields"
    ),
    {
        opts <- StepwiseOptions$new(
            model = "poisson",
            threshold = 0.05,
            direction = "both",
            frequency_col_name = "N_ID",
            capture_indicators = c("cap_1", "cap_2"),
            interaction_limit = 3
        )

        # Check inherited fields
        expect_equal(opts$model, "poisson")
        expect_equal(opts$threshold, 0.05)

        # Check StepwiseOptions-specific fields
        expect_equal(opts$direction, "both")
        expect_equal(opts$frequency_col_name, "N_ID")
        expect_equal(opts$capture_indicators, c("cap_1", "cap_2"))
        expect_equal(opts$interaction_limit, 3)
    }
)

# Test 4: EstimatorOptions constructor mapping
test_that(
    paste(
        "EstimatorOptions initializes all fields including",
        "inherited model and threshold"
    ),
    {
        opts <- EstimatorOptions$new(
            model = "logit",
            threshold = 0.10,
            list_pair = c("list_A", "list_B"),
            nfolds = 5,
            estimator = "ranger"
        )

        # Check inherited fields
        expect_equal(opts$model, "logit")
        expect_equal(opts$threshold, 0.10)

        # Check EstimatorOptions-specific fields
        expect_equal(opts$list_pair, c("list_A", "list_B"))
        expect_equal(opts$nfolds, 5)
        expect_equal(opts$estimator, "ranger")
    }
)

# Test 5: Class inheritance checks
test_that(
    paste("StepwiseOptions inherits from FrequencyOptions", "and Options"),
    {
        opts <- StepwiseOptions$new(
            model = "poisson",
            threshold = 0.05,
            direction = "forward",
            frequency_col_name = "N_ID",
            capture_indicators = c("cap_1", "cap_2")
        )

        expect_s3_class(opts, "StepwiseOptions")
        expect_s3_class(opts, "FrequencyOptions")
        expect_s3_class(opts, "Options")
    }
)

test_that(
    paste("AICOptions inherits from FrequencyOptions", "and Options"),
    {
        opts <- AICOptions$new(
            model = "poisson",
            frequency_col_name = "N_ID",
            capture_indicators = c("cap_1", "cap_2")
        )

        expect_s3_class(opts, "AICOptions")
        expect_s3_class(opts, "FrequencyOptions")
        expect_s3_class(opts, "Options")
    }
)

test_that(
    paste("EstimatorOptions inherits from Options"),
    {
        opts <- EstimatorOptions$new(
            model = "logit",
            threshold = 0.10,
            list_pair = c("list_A", "list_B"),
            nfolds = 5,
            estimator = "ranger"
        )

        expect_s3_class(opts, "EstimatorOptions")
        expect_s3_class(opts, "Options")
    }
)

# Test 6: Return-self behavior checks
test_that(
    paste("Options initialize method returns the same instance", "(self)"),
    {
        opts <- Options$new(model = "poisson", threshold = 0.05)
        result <- opts$initialize(model = "negbin", threshold = 0.10)

        expect_identical(result, opts)
        expect_equal(opts$model, "negbin")
        expect_equal(opts$threshold, 0.10)
    }
)

test_that(
    paste(
        "FrequencyOptions initialize, add_formula, add_formulas",
        "return the same instance"
    ),
    {
        opts <- FrequencyOptions$new(
            model = "poisson",
            threshold = 0.05,
            formulas = NULL,
            frequency_col_name = "N_ID"
        )
        init_result <- opts$initialize(
            model = "poisson",
            threshold = 0.05,
            formulas = NULL,
            frequency_col_name = "N_ID"
        )
        expect_identical(init_result, opts)

        formula1 <- stats::as.formula("N_ID ~ capture_1")
        add_result <- opts$add_formula(formula1)
        expect_identical(add_result, opts)

        formula2 <- stats::as.formula("N_ID ~ capture_2")
        add_multi_result <- opts$add_formulas(formula2)
        expect_identical(add_multi_result, opts)
    }
)

test_that(
    paste("StepwiseOptions initialize method returns the same", "instance"),
    {
        opts <- StepwiseOptions$new(
            model = "poisson",
            threshold = 0.05,
            direction = "forward",
            frequency_col_name = "N_ID",
            capture_indicators = c("cap_1", "cap_2")
        )
        result <- opts$initialize(
            model = "negbin",
            threshold = 0.10,
            direction = "backward",
            frequency_col_name = "count",
            capture_indicators = c("c1", "c2")
        )

        expect_identical(result, opts)
        expect_equal(opts$model, "negbin")
        expect_equal(opts$threshold, 0.10)
        expect_equal(opts$direction, "backward")
    }
)

test_that(
    paste("EstimatorOptions initialize method returns the same", "instance"),
    {
        opts <- EstimatorOptions$new(
            model = "logit",
            threshold = 0.10,
            list_pair = c("list_A", "list_B"),
            nfolds = 5,
            estimator = "ranger"
        )
        result <- opts$initialize(
            model = "gam",
            threshold = 0.15,
            list_pair = c("L1", "L2"),
            nfolds = 10,
            estimator = "logit"
        )

        expect_identical(result, opts)
        expect_equal(opts$model, "gam")
        expect_equal(opts$threshold, 0.15)
        expect_equal(opts$list_pair, c("L1", "L2"))
        expect_equal(opts$nfolds, 10)
    }
)
