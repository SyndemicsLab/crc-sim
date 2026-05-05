################################################################################
# File: test_crc_pipeline.R                                                    #
# Project: crc-sim                                                             #
# Created Date: 2026-05-05                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-05                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

test_that(
    paste(
        "create_data builds a binary contingency",
        "table with frequency counts"
    ),
    {
        set.seed(2026)
        data_table <- create_data(
            n_individuals = 500,
            p_captures = c(0.4, 0.3, 0.2),
            p_strata = 1
        )

        expected_columns <- c(
            "capture_1",
            "capture_2",
            "capture_3",
            "strata",
            "N_ID"
        )

        expect_s3_class(data_table, "data.table")
        expect_true(all(expected_columns %in% names(data_table)))
        expect_equal(sum(data_table$N_ID), 500)
        expect_true(all(data_table$N_ID >= 1))

        capture_columns <- c("capture_1", "capture_2", "capture_3")
        expect_true(all(unlist(data_table[, ..capture_columns]) %in% c(0, 1)))
    }
)

test_that("AIC model fitting supports poisson and negbin", {
    fixture_table <- data.table::data.table(
        capture_1 = c(1, 0, 1),
        capture_2 = c(0, 1, 1),
        N_ID = c(28, 24, 17)
    )

    required_columns <- c(
        "formula",
        "estimate",
        "AIC",
        "lower_ci",
        "upper_ci",
        "error"
    )

    poisson_result <- crc(
        data = fixture_table,
        freq_column = "N_ID",
        binary_variables = c("capture_1", "capture_2"),
        method = "poisson",
        formula_selection = "aic"
    )

    negbin_result <- suppressWarnings(crc(
        data = fixture_table,
        freq_column = "N_ID",
        binary_variables = c("capture_1", "capture_2"),
        method = "negbin",
        formula_selection = "aic"
    ))

    expect_true(all(required_columns %in% names(poisson_result)))
    expect_true(all(required_columns %in% names(negbin_result)))
    expect_true(anyNA(poisson_result$error))
    expect_true(anyNA(negbin_result$error))
    expect_true(any(is.na(poisson_result$error) & poisson_result$estimate > 0))
    expect_true(any(is.na(negbin_result$error) & negbin_result$estimate > 0))
})

test_that("stepwise model fitting supports poisson and negbin", {
    set.seed(2026)

    full_table <- create_data(
        n_individuals = 20000,
        p_captures = c(0.45, 0.35, 0.25, 0.2),
        p_strata = 1
    )

    capture_columns <- c(
        "capture_1",
        "capture_2",
        "capture_3",
        "capture_4"
    )

    observed_table <- full_table[rowSums(full_table[, ..capture_columns]) != 0]

    poisson_stepwise <- crc(
        data = observed_table,
        freq_column = "N_ID",
        binary_variables = capture_columns,
        method = "poisson",
        formula_selection = "stepwise",
        opts_stepwise = list(
            direction = "both",
            threshold = 0.05,
            verbose = FALSE
        )
    )

    negbin_stepwise <- suppressWarnings(crc(
        data = observed_table,
        freq_column = "N_ID",
        binary_variables = capture_columns,
        method = "negbin",
        formula_selection = "stepwise",
        opts_stepwise = list(
            direction = "both",
            threshold = 0.05,
            verbose = FALSE
        )
    ))

    expect_type(poisson_stepwise$estimate, "double")
    expect_type(negbin_stepwise$estimate, "double")
    expect_gt(poisson_stepwise$estimate, 0)
    expect_gt(negbin_stepwise$estimate, 0)
    expect_false(is.null(poisson_stepwise$formula))
    expect_false(is.null(negbin_stepwise$formula))
})

test_that(
    paste(
        "CRC estimates hidden zero-capture population",
        "with reasonable error"
    ),
    {
        set.seed(2026)

        full_table <- create_data(
            n_individuals = 20000,
            p_captures = c(0.45, 0.35, 0.25, 0.2),
            p_strata = 1
        )

        capture_columns <- c(
            "capture_1",
            "capture_2",
            "capture_3",
            "capture_4"
        )

        is_hidden <- rowSums(full_table[, ..capture_columns]) == 0
        ground_truth_hidden <- full_table[is_hidden, N_ID]
        observed_table <- full_table[!is_hidden]

        poisson_stepwise <- crc(
            data = observed_table,
            freq_column = "N_ID",
            binary_variables = capture_columns,
            method = "poisson",
            formula_selection = "stepwise",
            opts_stepwise = list(
                direction = "both",
                threshold = 0.05,
                verbose = FALSE
            )
        )

        negbin_stepwise <- suppressWarnings(crc(
            data = observed_table,
            freq_column = "N_ID",
            binary_variables = capture_columns,
            method = "negbin",
            formula_selection = "stepwise",
            opts_stepwise = list(
                direction = "both",
                threshold = 0.05,
                verbose = FALSE
            )
        ))

        poisson_relative_error <- abs(
            poisson_stepwise$estimate - ground_truth_hidden
        ) /
            ground_truth_hidden
        negbin_relative_error <- abs(
            negbin_stepwise$estimate - ground_truth_hidden
        ) /
            ground_truth_hidden

        expect_lt(poisson_relative_error, 0.35)
        expect_lt(negbin_relative_error, 0.35)
    }
)

test_that("run_crc output includes convergence timing benchmark column", {
    set.seed(2026)

    config <- list(
        fb0_05 = list(direction = "both", threshold = 0.05)
    )

    single_result <- run_crc(
        mode = "single",
        p_captures = c(0.45, 0.35, 0.25),
        p_strata = 1,
        n_individuals = 5000,
        methods = c("poisson"),
        config = config,
        suppress = 10,
        formula_selection = "stepwise",
        seed = 2026
    )

    bootstrap_result <- run_crc(
        mode = "bootstrap",
        n_bootstraps = 3,
        p_captures = c(0.45, 0.35, 0.25),
        p_strata = 1,
        n_individuals = 5000,
        methods = c("poisson"),
        config = config,
        suppress = 10,
        formula_selection = "stepwise",
        seed = 2026
    )

    expect_true("time_to_convergence_sec" %in% names(single_result))
    expect_true("time_to_convergence_sec" %in% names(bootstrap_result))
    expect_true(all(single_result$time_to_convergence_sec >= 0, na.rm = TRUE))
    expect_true(all(bootstrap_result$time_to_convergence_sec >= 0, na.rm = TRUE))
})
