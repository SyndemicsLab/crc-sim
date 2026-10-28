test_that("LoglinearOptions creates default AIC formulas", {
    opts <- LoglinearOptions$new(
        capture_columns = c("capture_1", "capture_2"),
        frequency_col_name = "N_ID"
    )

    expect_equal(opts$frequency_col_name, "N_ID")
    expect_equal(opts$model_family, "poisson")
    expect_equal(opts$selection_method, "aic")
    expect_equal(opts$selection_criterion, "AIC")
    expect_equal(
        opts$selection_options$formulas,
        formula_list("N_ID", c("capture_1", "capture_2"))
    )
})

test_that("LoglinearOptions stores stepwise settings with defaults", {
    opts <- LoglinearOptions$new(
        capture_columns = c("capture_1", "capture_2"),
        model_family = "negbin",
        selection_method = "stepwise",
        selection_criterion = "BIC",
        selection_options = list(direction = "forward")
    )

    expect_equal(opts$model_family, "negbin")
    expect_equal(opts$selection_method, "stepwise")
    expect_equal(opts$selection_criterion, "BIC")
    expect_equal(opts$selection_options$direction, "forward")
    expect_equal(opts$selection_options$interaction_limit, 2L)
})

test_that("LoglinearOptions rejects options for another selection method", {
    expect_error(
        LoglinearOptions$new(
            capture_columns = c("capture_1", "capture_2"),
            selection_options = list(direction = "both")
        ),
        "Unsupported selection option"
    )
})

test_that("LoglinearOptions rejects missing capture columns for generated formulas", {
    expect_error(
        LoglinearOptions$new(capture_columns = NULL),
        "capture_columns is required"
    )
})

test_that("LoglinearOptions rejects invalid model families", {
    expect_error(
        LoglinearOptions$new(
            capture_columns = c("capture_1", "capture_2"),
            model_family = "binomial"
        ),
        "should be one of"
    )
})

test_that("EstimatorOptions preserves estimator configuration", {
    opts <- EstimatorOptions$new(
        method = "plugin",
        capture_columns = c("capture_1", "capture_2"),
        threshold = 0.05,
        nuisance_function = "logit",
        nfolds = 5
    )

    expect_s3_class(opts, "EstimatorOptions")
    expect_s3_class(opts, "Options")
    expect_equal(opts$model, "plugin")
    expect_equal(opts$threshold, 0.05)
    expect_equal(opts$nuisance_function, "logit")
    expect_equal(opts$nfolds, 5)
})
