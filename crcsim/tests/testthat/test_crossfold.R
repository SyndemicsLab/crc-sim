################################################################################
# File: test_crossfold.R                                                       #
# Project: crcsim                                                              #
# Created Date: 2026-08-31                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("build_folds keeps all data in a single fold", {
    data <- data.frame(id = 1:4, value = letters[1:4])

    result <- build_folds(
        fold = 1,
        data = data,
        nfolds = 1,
        permutation = c(4, 3, 2, 1)
    )

    expect_named(result, c("train", "test", "sbset"))
    expect_equal(result$train, data)
    expect_equal(result$test, data)
    expect_identical(result$sbset, 1:4)
})

test_that("build_folds partitions permuted rows into balanced folds", {
    data <- data.frame(id = 1:6)
    permutation <- c(6, 1, 4, 2, 5, 3)

    first <- build_folds(1, data, nfolds = 3, permutation = permutation)
    second <- build_folds(2, data, nfolds = 3, permutation = permutation)
    third <- build_folds(3, data, nfolds = 3, permutation = permutation)

    expect_identical(first$sbset, 1:2)
    expect_identical(second$sbset, 3:4)
    expect_identical(third$sbset, 5:6)
    expect_equal(first$test, c(6, 1))
    expect_equal(second$test, c(4, 2))
    expect_equal(third$test, c(5, 3))
    expect_equal(sort(c(first$test, second$test, third$test)), 1:6)
})

test_that("build_folds puts remainder rows in the final fold", {
    data <- data.frame(id = 1:5)
    permutation <- c(5, 4, 3, 2, 1)

    first <- build_folds(1, data, nfolds = 2, permutation = permutation)
    second <- build_folds(2, data, nfolds = 2, permutation = permutation)

    expect_length(first$test, 3)
    expect_length(second$test, 2)
    expect_equal(first$test, c(5, 4, 3))
    expect_equal(second$test, c(2, 1))
    expect_equal(sort(c(first$test, second$test)), 1:5)
})

test_that("run_fold returns NULL when nuisance estimation fails", {
    fold_split <- list(
        train = data.frame(
            capture_1 = c(0, 1, 0, 1),
            capture_2 = c(0, 0, 1, 1)
        ),
        test = data.frame(capture_1 = c(0, 1), capture_2 = c(1, 1)),
        sbset = 1:2
    )

    result <- run_fold(
        fold_split = fold_split,
        method = "plugin",
        j = 1,
        k = 2,
        margin = 0.05,
        n_lists = 99,
        n = 4,
        func = "logit"
    )

    expect_null(result)
})
