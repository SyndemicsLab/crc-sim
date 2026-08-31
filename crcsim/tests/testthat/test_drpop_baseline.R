################################################################################
# File: test_drpop_baseline.R                                                  #
# Project: crcsim                                                              #
# Created Date: 2026-08-31                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-31                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

make_drpop_baseline_data <- function(
    n_individuals = 300,
    n_covariates = 1,
    captures = paste0("capture_", seq_len(2))
) {
    set.seed(20260831)
    covariate_columns <- if (n_covariates == 0) {
        character()
    } else {
        paste0("covariate_", seq_len(n_covariates))
    }
    covariate_function <- function(data_table) {
        if (n_covariates == 0) {
            return(rep(0.5, nrow(data_table)))
        }
        return(stats::plogis(data_table$covariate_1))
    }
    baseline_data <- build_covariate_data(n_individuals, n_covariates) |>
        captures_from_covariates(
            length(captures),
            covariate_function
        ) |>
        relocate_covariates(covariate_columns) |>
        all_int_cols_to_numeric() |>
        extract_captured_data(captures) |>
        as.data.frame()
    return(baseline_data)
}

run_drpop_baseline <- function(
    data,
    captures = paste0("capture_", seq_len(2)),
    nfolds = 1,
    margin = 0.005
) {
    nuisance_fit <- drpop::popsize(
        data = data,
        K = length(captures),
        capturelists = as.numeric(seq_along(captures)),
        funcname = "logit",
        nfolds = nfolds,
        margin = margin
    )
    if (length(captures) == 2 && ncol(data) == length(captures)) {
        return(nuisance_fit$result)
    }
    result <- drpop::popsize(
        data = data,
        capturelists = as.numeric(seq_along(captures)),
        getnuis = nuisance_fit$nuis,
        idfold = nuisance_fit$idfold
    )$result
    return(result)
}

expect_crcsim_matches_drpop <- function(data, method, drpop_method) {
    drpop_result <- suppressWarnings(run_drpop_baseline(
        data,
        captures = c("capture_1", "capture_2"),
        nfolds = 1,
        margin = 0.005
    ))
    crcsim_result <- suppressWarnings(estimate_capture_prob(
        data,
        n_lists = 2,
        method = method,
        func = "logit",
        nfolds = 1,
        margin = 0.005,
        seed = 20260831
    ))

    expect_true(all(c("method", "n", "cin.l", "cin.u") %in%
        names(drpop_result)))
    drpop_row <- match(drpop_method, toupper(drpop_result$method))
    expect_false(
        is.na(drpop_row),
        info = paste("Available drpop methods:",
            toString(drpop_result$method))
    )

    expect_equal(
        crcsim_result$n,
        drpop_result$n[[drpop_row]],
        tolerance = 1
    )
    expect_equal(
        crcsim_result$ci_l,
        drpop_result$cin.l[[drpop_row]],
        tolerance = 1
    )
    expect_equal(
        crcsim_result$ci_u,
        drpop_result$cin.u[[drpop_row]],
        tolerance = 1
    )
    return(invisible(NULL))
}

test_that("drpop and crcsim agree without covariates", {
    skip_if_not_installed("drpop")
    data <- make_drpop_baseline_data(n_covariates = 0)

    drpop_result <- suppressWarnings(run_drpop_baseline(
        data,
        captures = c("capture_1", "capture_2"),
        nfolds = 1
    ))
    crcsim_result <- estimate_capture_prob(
        data,
        n_lists = 2,
        method = "plugin",
        func = "logit",
        nfolds = 1,
        margin = 0.005
    )

    expect_true("n" %in% names(drpop_result))
    expect_named(
        crcsim_result,
        c("listpair", "n", "sigma_n", "ci_l", "ci_u")
    )
    expect_equal(crcsim_result$listpair, "1,2")
    expect_true(any(abs(drpop_result$n - crcsim_result$n) <= 1))
})

test_that("drpop and crcsim return comparable covariate estimates", {
    skip_if_not_installed("drpop")
    data <- make_drpop_baseline_data()

    drpop_result <- suppressWarnings(run_drpop_baseline(
        data,
        captures = c("capture_1", "capture_2"),
        nfolds = 1,
        margin = 0.005
    ))
    crcsim_result <- suppressWarnings(estimate_capture_prob(
        data,
        n_lists = 2,
        method = "plugin",
        func = "logit",
        nfolds = 1,
        margin = 0.005,
        seed = 20260831
    ))

    expect_true(all(c("method", "n", "cin.l", "cin.u") %in%
        names(drpop_result)))
    expect_named(
        crcsim_result,
        c("listpair", "n", "sigma_n", "ci_l", "ci_u")
    )
    expect_equal(crcsim_result$listpair, "1,2")
    expect_true(all(is.finite(unlist(crcsim_result[1, -1]))))
    expect_true(all(is.finite(drpop_result$n)))
})

test_that("drpop and crcsim return comparable doubly robust estimates", {
    skip_if_not_installed("drpop")
    data <- make_drpop_baseline_data()

    drpop_result <- suppressWarnings(run_drpop_baseline(
        data,
        captures = c("capture_1", "capture_2"),
        nfolds = 1,
        margin = 0.005
    ))
    crcsim_result <- suppressWarnings(estimate_capture_prob(
        data,
        n_lists = 2,
        method = "doubly_robust",
        func = "logit",
        nfolds = 1,
        margin = 0.005,
        seed = 20260831
    ))

    expect_named(
        crcsim_result,
        c("listpair", "n", "sigma_n", "ci_l", "ci_u")
    )
    expect_equal(crcsim_result$listpair, "1,2")
    expect_true(all(is.finite(unlist(crcsim_result[1, -1]))))
    expect_equal(crcsim_result$n, drpop_result$n[[2]], tolerance = 1)
})

test_that("drpop and crcsim return comparable TMLE estimates", {
    skip_if_not_installed("drpop")
    data <- make_drpop_baseline_data()

    drpop_result <- suppressWarnings(run_drpop_baseline(
        data,
        captures = c("capture_1", "capture_2"),
        nfolds = 1,
        margin = 0.005
    ))
    crcsim_result <- suppressWarnings(estimate_capture_prob(
        data,
        n_lists = 2,
        method = "tmle",
        func = "logit",
        nfolds = 1,
        margin = 0.005,
        seed = 20260831
    ))

    expect_named(
        crcsim_result,
        c("listpair", "n", "sigma_n", "ci_l", "ci_u")
    )
    expect_equal(crcsim_result$listpair, "1,2")
    expect_true(all(is.finite(unlist(crcsim_result[1, -1]))))
    expect_true(any(abs(drpop_result$n - crcsim_result$n) <= 1))
})

test_that("crcsim plugin matches drpop PLUGIN", {
    skip_if_not_installed("drpop")
    expect_crcsim_matches_drpop(
        make_drpop_baseline_data(),
        method = "plugin",
        drpop_method = "PI"
    )
})

test_that("crcsim doubly robust matches drpop DR", {
    skip_if_not_installed("drpop")
    expect_crcsim_matches_drpop(
        make_drpop_baseline_data(),
        method = "doubly_robust",
        drpop_method = "DR"
    )
})

test_that("crcsim TMLE matches drpop TMLE", {
    skip_if_not_installed("drpop")
    expect_crcsim_matches_drpop(
        make_drpop_baseline_data(),
        method = "tmle",
        drpop_method = "TMLE"
    )
})
