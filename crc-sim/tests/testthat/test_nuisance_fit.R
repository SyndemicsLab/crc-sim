################################################################################
# File: test_nuisance_fit.R                                                    #
# Project: crc-sim                                                             #
# Created Date: 2026-05-20                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-20                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

# ============================================================================
# Test Fixtures and Helpers
# ============================================================================

#' Create row-level data with binary capture columns and covariates
make_nuisance_fixture <- function(
    n_individuals = 200,
    n_captures = 2,
    with_covariates = TRUE
) {
    set.seed(42)
    data <- as.data.frame(create_data(
        n_individuals = n_individuals,
        n_captures = n_captures,
        return_frequency_table = FALSE
    ))

    if (with_covariates) {
        # Add a continuous covariate and a binary covariate
        data$age <- rnorm(nrow(data), mean = 50, sd = 15)
        data$gender <- rbinom(nrow(data), 1, 0.5)
    }

    return(data)
}

#' Create train/test splits from row-level data
make_train_test_split <- function(data, train_fraction = 0.7) {
    set.seed(123)
    n_rows <- nrow(data)
    train_idx <- sample(1:n_rows, size = ceiling(n_rows * train_fraction))
    test_idx <- setdiff(1:n_rows, train_idx)

    return(list(
        train = data[train_idx, , drop = FALSE],
        test = data[test_idx, , drop = FALSE]
    ))
}

#' Create minimal data with no covariates for fallback path testing
make_minimal_fixture <- function() {
    set.seed(99)
    data <- data.frame(
        capture_1 = c(1, 0, 1, 1, 0),
        capture_2 = c(1, 1, 0, 1, 1)
    )
    return(data)
}

# ============================================================================
# empirical_means Tests
# ============================================================================

test_that("empirical_means returns list with expected names", {
    data <- make_minimal_fixture()
    train <- data[1:3, ]
    test <- data[4:5, ]

    result <- crcsim:::empirical_means(train, test, j = 1, k = 2)

    expect_type(result, "list")
    expect_named(result, c("q_j", "q_k", "q_jk"))
})

test_that("empirical_means returns vectors matching test row count", {
    data <- make_minimal_fixture()
    train <- data[1:3, ]
    test <- data[4:5, ]

    result <- crcsim:::empirical_means(train, test, j = 1, k = 2)

    expect_equal(length(result$q_j), nrow(test))
    expect_equal(length(result$q_k), nrow(test))
    expect_equal(length(result$q_jk), nrow(test))
})

test_that("empirical_means values are constant repeats of training means", {
    data <- make_minimal_fixture()
    train <- data[1:3, ]
    test <- data[4:5, ]

    result <- crcsim:::empirical_means(train, test, j = 1, k = 2)

    expected_q_j <- mean(train[[1]])
    expected_q_k <- mean(train[[2]])
    expected_q_jk <- mean(train[[1]] * train[[2]])

    expect_true(all(result$q_j == expected_q_j))
    expect_true(all(result$q_k == expected_q_k))
    expect_true(all(result$q_jk == expected_q_jk))
})

test_that("empirical_means handles minimal data sizes", {
    # Single row train/test
    train <- data.frame(capture_1 = 1, capture_2 = 0)
    test <- data.frame(capture_1 = 1, capture_2 = 1)

    result <- crcsim:::empirical_means(train, test, j = 1, k = 2)

    expect_equal(nrow(test), length(result$q_j))
    expect_equal(nrow(test), length(result$q_k))
    expect_equal(nrow(test), length(result$q_jk))
})

# ============================================================================
# fit_glm Tests
# ============================================================================

test_that("fit_glm returns object inheriting glm", {
    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    train <- splits$train

    formula_obj <- stats::formula("capture_1 ~ age + gender")
    model <- crcsim:::fit_glm(train, formula_obj)

    expect_s3_class(model, "glm")
})

test_that("fit_glm binomial link is logit", {
    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    train <- splits$train

    formula_obj <- stats::formula("capture_1 ~ age + gender")
    model <- crcsim:::fit_glm(train, formula_obj)

    expect_equal(model$family$link, "logit")
})

# ============================================================================
# fit_gam Tests
# ============================================================================

test_that("fit_gam returns object inheriting gam", {
    skip_if_not_installed("mgcv")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    train <- splits$train

    formula_obj <- stats::formula("capture_1 ~ age + gender")
    model <- crcsim:::fit_gam(train, formula_obj)

    expect_s3_class(model, "gam")
})

# ============================================================================
# fit_random_forest Tests
# ============================================================================

test_that("fit_random_forest returns ranger object", {
    skip_if_not_installed("ranger")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    train <- splits$train

    formula_obj <- stats::formula("factor(capture_1) ~ age + gender")
    model <- crcsim:::fit_random_forest(train, formula_obj)

    expect_s3_class(model, "ranger")
})

test_that("fit_random_forest uses probability classification", {
    skip_if_not_installed("ranger")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    train <- splits$train

    formula_obj <- stats::formula("factor(capture_1) ~ age + gender")
    model <- crcsim:::fit_random_forest(train, formula_obj)

    expect_true(model$probability)
})

# ============================================================================
# predict_qhat Tests
# ============================================================================

test_that("predict_qhat returns numeric vector of test row length", {
    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    train <- splits$train
    test <- splits$test

    formula_obj <- stats::formula("capture_1 ~ age + gender")
    model <- crcsim:::fit_glm(train, formula_obj)
    predictions <- crcsim:::predict_qhat(model, test)

    expect_type(predictions, "double")
    expect_equal(length(predictions), nrow(test))
})

test_that("predict_qhat predictions are between 0 and 1 for glm", {
    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    train <- splits$train
    test <- splits$test

    formula_obj <- stats::formula("capture_1 ~ age + gender")
    model <- crcsim:::fit_glm(train, formula_obj)
    predictions <- crcsim:::predict_qhat(model, test)

    expect_true(all(predictions >= 0 & predictions <= 1))
})

# ============================================================================
# glm_qhats Tests
# ============================================================================

test_that("glm_qhats no-covariate branch returns fallback output", {
    # Create data with only capture columns (no covariates)
    train <- data.frame(capture_1 = c(1, 0, 1), capture_2 = c(1, 1, 0))
    test <- data.frame(capture_1 = c(0, 1), capture_2 = c(1, 0))

    result <- crcsim:::glm_qhats(
        train = train,
        test = test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_named(result, c("q_j", "q_k", "q_jk", "y_j", "y_k"))
})

test_that("glm_qhats covariate branch returns list with expected names", {
    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    result <- glm_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_named(result, c("q_j", "q_k", "q_jk", "y_j", "y_k"))
})

test_that("glm_qhats output vector lengths match test rows", {
    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    result <- glm_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_equal(length(result$q_j), nrow(splits$test))
    expect_equal(length(result$q_k), nrow(splits$test))
    expect_equal(length(result$q_jk), nrow(splits$test))
})

test_that("glm_qhats q_j/q_k/q_jk respect margin lower bound", {
    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    margin <- 0.01

    result <- glm_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = margin
    )

    expect_true(all(result$q_j >= margin))
    expect_true(all(result$q_k >= margin))
    expect_true(all(result$q_jk >= margin))
})

test_that("glm_qhats y_j and y_k match test capture columns exactly", {
    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    result <- glm_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_equal(result$y_j, splits$test[[1]])
    expect_equal(result$y_k, splits$test[[2]])
})

# ============================================================================
# gam_qhats Tests
# ============================================================================

test_that("gam_qhats no-covariate branch returns fallback output", {
    skip_if_not_installed("mgcv")

    # Create data with only capture columns (no covariates)
    train <- data.frame(capture_1 = c(1, 0, 1), capture_2 = c(1, 1, 0))
    test <- data.frame(capture_1 = c(0, 1), capture_2 = c(1, 0))

    result <- crcsim:::gam_qhats(
        train = train,
        test = test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_named(result, c("q_j", "q_k", "q_jk", "y_j", "y_k"))
})

test_that("gam_qhats covariate branch returns list with expected names", {
    skip_if_not_installed("mgcv")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    result <- gam_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_named(result, c("q_j", "q_k", "q_jk", "y_j", "y_k"))
})

test_that("gam_qhats output vector lengths match test rows", {
    skip_if_not_installed("mgcv")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    result <- gam_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_equal(length(result$q_j), nrow(splits$test))
    expect_equal(length(result$q_k), nrow(splits$test))
    expect_equal(length(result$q_jk), nrow(splits$test))
})

test_that("gam_qhats q_j/q_k/q_jk respect margin lower bound", {
    skip_if_not_installed("mgcv")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    margin <- 0.01

    result <- gam_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = margin
    )

    expect_true(all(result$q_j >= margin))
    expect_true(all(result$q_k >= margin))
    expect_true(all(result$q_jk >= margin))
})

test_that("gam_qhats y_j and y_k match test capture columns exactly", {
    skip_if_not_installed("mgcv")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    result <- gam_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_equal(result$y_j, splits$test[[1]])
    expect_equal(result$y_k, splits$test[[2]])
})

# ============================================================================
# random_forest_qhats Tests
# ============================================================================

test_that("random_forest_qhats no-covariate branch returns fallback output", {
    skip_if_not_installed("ranger")

    # Create data with only capture columns (no covariates)
    train <- data.frame(capture_1 = c(1, 0, 1), capture_2 = c(1, 1, 0))
    test <- data.frame(capture_1 = c(0, 1), capture_2 = c(1, 0))

    result <- crcsim:::random_forest_qhats(
        train = train,
        test = test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_named(result, c("q_j", "q_k", "q_jk", "y_j", "y_k"))
})

test_that("random_forest_qhats covariate branch returns list with expected names", {
    skip_if_not_installed("ranger")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    result <- random_forest_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_named(result, c("q_j", "q_k", "q_jk", "y_j", "y_k"))
})

test_that("random_forest_qhats output vector lengths match test rows", {
    skip_if_not_installed("ranger")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    result <- random_forest_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_equal(length(result$q_j), nrow(splits$test))
    expect_equal(length(result$q_k), nrow(splits$test))
    expect_equal(length(result$q_jk), nrow(splits$test))
})

test_that("random_forest_qhats q_j/q_k/q_jk respect margin lower bound", {
    skip_if_not_installed("ranger")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)
    margin <- 0.01

    result <- random_forest_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = margin
    )

    expect_true(all(result$q_j >= margin))
    expect_true(all(result$q_k >= margin))
    expect_true(all(result$q_jk >= margin))
})

test_that("random_forest_qhats y_j and y_k match test capture columns exactly", {
    skip_if_not_installed("ranger")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    result <- random_forest_qhats(
        train = splits$train,
        test = splits$test,
        n_capture_cols = 2,
        j = 1,
        k = 2,
        margin = 0.005
    )

    expect_equal(result$y_j, splits$test[[1]])
    expect_equal(result$y_k, splits$test[[2]])
})

# ============================================================================
# Error and Robustness Tests
# ============================================================================

test_that("glm_qhats errors on mismatched train/test columns", {
    train <- data.frame(
        capture_1 = c(1, 0),
        capture_2 = c(1, 1),
        age = c(50, 60)
    )
    test <- data.frame(capture_1 = c(1, 0), capture_2 = c(1, 1))
    # test is missing 'age' column

    expect_error(
        crcsim:::glm_qhats(
            train = train,
            test = test,
            n_capture_cols = 2,
            j = 1,
            k = 2
        ),
        "identical"
    )
})

test_that("gam_qhats errors on mismatched train/test columns", {
    skip_if_not_installed("mgcv")

    train <- data.frame(
        capture_1 = c(1, 0),
        capture_2 = c(1, 1),
        age = c(50, 60)
    )
    test <- data.frame(capture_1 = c(1, 0), capture_2 = c(1, 1))

    expect_error(
        crcsim:::gam_qhats(
            train = train,
            test = test,
            n_capture_cols = 2,
            j = 1,
            k = 2
        ),
        "identical"
    )
})

test_that("random_forest_qhats errors on mismatched train/test columns", {
    skip_if_not_installed("ranger")

    train <- data.frame(
        capture_1 = c(1, 0),
        capture_2 = c(1, 1),
        age = c(50, 60)
    )
    test <- data.frame(capture_1 = c(1, 0), capture_2 = c(1, 1))

    expect_error(
        crcsim:::random_forest_qhats(
            train = train,
            test = test,
            n_capture_cols = 2,
            j = 1,
            k = 2
        ),
        "identical"
    )
})

test_that("glm_qhats with invalid j index errors on column access", {
    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    expect_error(
        crcsim:::glm_qhats(
            train = splits$train,
            test = splits$test,
            n_capture_cols = 2,
            j = 99, # Invalid index
            k = 2,
            margin = 0.005
        )
    )
})

test_that("gam_qhats with invalid k index errors on column access", {
    skip_if_not_installed("mgcv")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    expect_error(
        crcsim:::gam_qhats(
            train = splits$train,
            test = splits$test,
            n_capture_cols = 2,
            j = 1,
            k = 99, # Invalid index
            margin = 0.005
        )
    )
})

test_that("random_forest_qhats with invalid j index errors on column access", {
    skip_if_not_installed("ranger")

    data <- make_nuisance_fixture(with_covariates = TRUE)
    splits <- make_train_test_split(data)

    expect_error(
        crcsim:::random_forest_qhats(
            train = splits$train,
            test = splits$test,
            n_capture_cols = 2,
            j = 99, # Invalid index
            k = 2,
            margin = 0.005
        )
    )
})
