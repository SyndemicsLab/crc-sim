################################################################################
# File: temp.R                                                                 #
# Project: crcsim                                                              #
# Created Date: 2026-08-19                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-24                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

devtools::load_all("crcsim")
library(dplyr)
library(purrr)

set.seed(2024)
n_individuals <- 30000
n_sims <- 100

captures <- c(
    "capture_1",
    "capture_2",
    "capture_3",
    "capture_4",
    "capture_5",
    "capture_6"
)

crc_loop <- function(intercept) {
    n_covariates <- 4

    # range from -2 to 2 (incr. by 0.4)
    # intercept <- 0.0

    # change regression coefficients to -1.5 and 1.5
    regression_coef <- 0

    log_regress_covariates <- function(data_table) {
        z1_values <- data_table[["covariate_1"]]
        z2_values <- data_table[["covariate_2"]]
        z3_values <- data_table[["covariate_3"]]
        z4_values <- data_table[["covariate_4"]]

        v <- exp(
            intercept +
                regression_coef * z1_values +
                regression_coef * z2_values -
                regression_coef * z3_values -
                regression_coef * z4_values
        )

        capture_prob <- v / (1 + v)
        return(capture_prob)
    }

    # nolint start: implicit_assignment_linter
    # Necessary to suppress output from the simulations. The results are captured
    # and stored in the `res` variable.
    invisible(capture.output(
        res <- map(
            seq_len(n_sims),
            ~ crcsim::simulate(
                n_individuals = n_individuals,
                n_covariates = n_covariates,
                captures = captures,
                cov_func = log_regress_covariates
            )
        )
    ))
    # nolint end: implicit_assignment_linter

    final_results <- bind_rows(res) |>
        mutate(coverage = as.numeric(lower_ci <= 0 & upper_ci >= 0)) |>
        group_by(method) |>
        summarize(
            estimate = mean(estimate),
            lower_ci = mean(lower_ci),
            upper_ci = mean(upper_ci),
            coverage = sum(coverage) / n()
        )

    name <- paste0("final_results_", intercept, "int_0.0rc.csv")
    write.csv(
        final_results,
        name,
        row.names = FALSE
    )
}

lapply(
    seq(-2, 2, by = 0.4),
    crc_loop
)
