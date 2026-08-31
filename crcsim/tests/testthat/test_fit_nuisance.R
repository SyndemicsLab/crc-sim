################################################################################
# File: test_fit_nuisance.R                                                    #
# Project: crcsim                                                              #
# Created Date: 2026-08-31                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("qhat_logit returns bounded nuisance estimates", {
    train <- data.frame(
        capture_1 = c(0, 1, 0, 1, 1, 0, 1, 0),
        capture_2 = c(0, 0, 1, 1, 1, 0, 0, 1),
        covariate = c(0, 1, 0, 1, 0, 1, 0, 1)
    )
    test <- train[1:3, ]

    result <- crcsim:::qhat_logit(
        train,
        test,
        n_lists = 2,
        j = 1,
        k = 2,
        margin = 0.05
    )

    expect_named(result, c("q_1", "q_2", "q_12"))
    expect_length(result$q_1, nrow(test))
    expect_length(result$q_2, nrow(test))
    expect_length(result$q_12, nrow(test))
    expect_true(all(result$q_1 >= 0.05 & result$q_1 <= 1))
    expect_true(all(result$q_2 >= 0.05 & result$q_2 <= 1))
    expect_true(all(result$q_12 >= 0.05 & result$q_12 <= 1))
})

test_that("qhat_logit handles repeated person-level rows", {
    person_data <- data.frame(
        capture_1 = c(0, 1, 0, 1, 1, 0, 1, 0),
        capture_2 = c(0, 0, 1, 1, 1, 0, 0, 1),
        covariate = c(0, 1, 0, 1, 0, 1, 0, 1)
    )
    expanded_data <- person_data[rep(seq_len(nrow(person_data)), 1000), ]

    result <- crcsim:::qhat_logit(
        expanded_data,
        expanded_data[1:3, ],
        n_lists = 2,
        j = 1,
        k = 2,
        margin = 0.05
    )

    expect_named(result, c("q_1", "q_2", "q_12"))
    expect_length(result$q_1, 3)
})

test_that("qhat_logit returns NULL when a GLM cannot be fitted", {
    train <- data.frame(
        capture_1 = c("no", "yes", "no"),
        capture_2 = c("no", "no", "yes"),
        covariate = c(0, 1, 0)
    )
    test <- train

    expect_warning(
        result <- crcsim:::qhat_logit(
            train,
            test,
            n_lists = 2,
            j = 1,
            k = 2,
            margin = 0.05
        ),
        "One or more GLM fits failed"
    )
    expect_null(result)
})

test_that("nuisance_estimation delegates to qhat_logit", {
    train <- data.frame(
        capture_1 = c(0, 1, 0, 1),
        capture_2 = c(0, 0, 1, 1),
        covariate = c(0, 1, 0, 1)
    )
    test <- train

    direct <- crcsim:::qhat_logit(train, test, 2, 1, 2, 0.05)
    delegated <- crcsim:::nuisance_estimation(
        "logit",
        train,
        test,
        2,
        1,
        2,
        0.05
    )

    expect_identical(delegated, direct)
})

test_that("tmle_nuisance preserves initial values with zero iterations", {
    result <- crcsim:::tmle_nuisance(
        q_1 = c(0.6, 0.7),
        q_2 = c(0.7, 0.8),
        q_12 = c(0.2, 0.3),
        y_j = c(0, 1),
        y_k = c(1, 1),
        y_jk = c(0, 1),
        iterations = 0,
        margin = 0.1,
        n_lists = 3
    )

    expect_identical(names(result), c("q_1", "q_2", "q_12"))
    expect_equal(result$q_1, c(0.6, 0.7))
    expect_equal(result$q_2, c(0.7, 0.8))
    expect_equal(result$q_12, c(0.2, 0.3))
})

test_that("tmle_nuisance uses the two-list complement for q_2", {
    result <- crcsim:::tmle_nuisance(
        q_1 = 0.6,
        q_2 = 0.7,
        q_12 = 0.2,
        y_j = 0,
        y_k = 1,
        y_jk = 0,
        iterations = 1,
        margin = 0.1,
        n_lists = 2
    )

    expect_equal(result$q_2 + result$q_1 - result$q_12, 1)
})

test_that("fit_glm_tmle handles complete cases", {
    result <- crcsim:::fit_glm_tmle(
        response = c(0, 1, NA, 1),
        offset = c(0, 0, 0, NA),
        ratio = c(1, 2, 3, 4)
    )

    expect_named(result, c("value", "error"))
    expect_length(result$value, 2)
    expect_true(all(is.finite(result$value)))
    expect_length(result$error, 1)
    expect_true(is.finite(result$error))
})

test_that("fit_glm_tmle reports failed fits", {
    result <- crcsim:::fit_glm_tmle(
        response = c(NA, NA),
        offset = c(0, 0),
        ratio = c(1, 1)
    )

    expect_null(result$value)
    expect_identical(result$error, Inf)
})
