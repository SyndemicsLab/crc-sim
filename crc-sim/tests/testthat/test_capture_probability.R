################################################################################
# File: test_capture_probability.R                                             #
# Project: crc-sim                                                             #
# Created Date: 2026-05-22                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-26                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that("scalar conditional_capture_prob returns correct gamma_hat", {
    q_j <- 0.5
    q_k <- 0.4
    q_jk <- 0.2
    expected_gamma_hat <- q_jk / (q_j * q_k)
    gamma_hat <- conditional_capture_prob(q_j, q_k, q_jk)
    expect_equal(gamma_hat, expected_gamma_hat)
})

test_that("vectorized conditional_capture_prob returns correct gamma_hat", {
    q_j <- c(0.5, 0.3, 0.6)
    q_k <- c(0.4, 0.2, 0.5)
    q_jk <- c(0.2, 0.1, 0.3)
    expected_gamma_hat <- q_jk / (q_j * q_k)
    gamma_hat <- conditional_capture_prob(q_j, q_k, q_jk)
    expect_equal(gamma_hat, expected_gamma_hat)
})

test_that("scalar plugin_estimation returns correct psi_hat", {
    q_j <- 0.5
    q_k <- 0.4
    q_jk <- 0.2
    gamma_hat <- q_jk / (q_j * q_k)
    expected_psi_hat <- 1 / mean(1 / gamma_hat)
    psi_hat <- get_plugin_estimation(q_j, q_k, q_jk)
    expect_equal(psi_hat, expected_psi_hat)
})

test_that("vectorized plugin_estimation returns correct psi_hat", {
    q_j <- c(0.5, 0.3, 0.6)
    q_k <- c(0.4, 0.2, 0.5)
    q_jk <- c(0.2, 0.1, 0.3)
    gamma_hat <- q_jk / (q_j * q_k)
    expected_psi_hat <- 1 / mean(1 / gamma_hat)
    psi_hat <- get_plugin_estimation(q_j, q_k, q_jk)
    expect_equal(psi_hat, expected_psi_hat)
})
