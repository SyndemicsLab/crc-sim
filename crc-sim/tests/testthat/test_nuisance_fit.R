################################################################################
# File: test_nuisance_fit.R                                                    #
# Project: crc-sim                                                             #
# Created Date: 2026-05-20                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-22                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("verify empirical_qhats is correct", {
    set.seed(123)
    sim_data <- create_data(
        n_individuals = 300,
        n_captures = 3
    )
    train <- sim_data[1:200, ]
    test <- sim_data[201:300, ]
    qhats <- empirical_qhats(train, test, 1, 2)

    q_j <- rep(mean(train[[1]]), nrow(test))
    q_k <- rep(mean(train[[2]]), nrow(test))
    q_jk <- rep(mean(train[[1]] * train[[2]]), nrow(test))
    y_j <- test[[1]]
    y_k <- test[[2]]

    expect_equal(qhats$q_j, q_j)
    expect_equal(qhats$q_k, q_k)
    expect_equal(qhats$q_jk, q_jk)
    expect_equal(qhats$y_j, y_j)
    expect_equal(qhats$y_k, y_k)
})

test_that("verify empirical_qhats handles invalid j values", {
    set.seed(123)
    sim_data <- create_data(
        n_individuals = 300,
        n_captures = 3
    )
    train <- sim_data[1:200, ]
    test <- sim_data[201:300, ]
    expect_error(empirical_qhats(train, test, 5, 2))
})

test_that("verify empirical_qhats handles invalid k values", {
    set.seed(123)
    sim_data <- create_data(
        n_individuals = 300,
        n_captures = 3
    )
    train <- sim_data[1:200, ]
    test <- sim_data[201:300, ]
    expect_error(empirical_qhats(train, test, 1, 5))
})

test_that("verify we get qhats from glm_qhats", {
    set.seed(1234)
    sim_data <- create_data(
        n_individuals = 300,
        n_captures = 3
    )
    train <- sim_data[1:200, ]
    test <- sim_data[201:300, ]
    qhats <- glm_qhats(train, test, 3, 1, 2)

    q_j <- rep(0.255, nrow(test))
    names(q_j) <- seq_len(nrow(test))
    q_k <- rep(0.31, nrow(test))
    names(q_k) <- seq_len(nrow(test))
    q_jk <- rep(0.22, nrow(test))
    names(q_jk) <- seq_len(nrow(test))

    expect_equal(qhats$q_j, q_j)
    expect_equal(qhats$q_k, q_k)
    expect_equal(qhats$q_jk, q_jk)
    expect_equal(qhats$y_j, test[[1]])
    expect_equal(qhats$y_k, test[[2]])
})

test_that("verify we get qhats from gam_qhats", {
    set.seed(1234)
    sim_data <- create_data(
        n_individuals = 300,
        n_captures = 3
    )
    train <- sim_data[1:200, ]
    test <- sim_data[201:300, ]
    qhats <- gam_qhats(train, test, 3, 1, 2)

    q_j <- array(rep(0.475, nrow(test)))
    q_k <- array(rep(0.53, nrow(test)))
    q_jk <- array(rep(0.22, nrow(test)))

    expect_equal(unname(qhats$q_j), q_j)
    expect_equal(unname(qhats$q_k), q_k)
    expect_equal(unname(qhats$q_jk), q_jk)
    expect_equal(qhats$y_j, test[[1]])
    expect_equal(qhats$y_k, test[[2]])
})
