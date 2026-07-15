################################################################################
# File: sim.R                                                                  #
# Project: crc-sim                                                             #
# Created Date: 2026-06-18                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-07-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

library(crcsim)
library(drpop)
library(dplyr)
library(purrr)

run_crc_simulations <- function(
    n_individuals,
    n_captures,
    p_captures,
    covariate_df,
    captures
) {
    sim_data <- create_data(
        n_individuals = n_individuals,
        n_captures = n_captures,
        p_captures = p_captures,
        covariate_ranges = covariate_df
    ) |>
        all_int_cols_to_numeric() |>
        extract_captured_data(captures)

    no_cov_sim_data <- sim_data |>
        select(-covariate_1, -covariate_2) |>
        build_contingency_table()

    known_pop_size <- nrow(sim_data)
    unknown_pop_size <- n_individuals - known_pop_size

    qhat <- drpop::popsize(
        data = sim_data,
        funcname = c("logit"),
        nfolds = 2,
        margin = 0.005
    )

    psin_estimates <- drpop::popsize(
        data = sim_data,
        getnuis = qhat$nuis,
        idfold = qhat$idfold
    )$result

    drpop_df <- dplyr::select(
        psin_estimates,
        c("method", "n", "cin.l", "cin.u")
    ) |>
        dplyr::rename(
            estimate = n,
            lower_ci = cin.l,
            upper_ci = cin.u
        ) |>
        dplyr::mutate(
            estimate = ((estimate - n_individuals) / n_individuals) * 100,
            lower_ci = ((lower_ci - n_individuals) / n_individuals) * 100,
            upper_ci = ((upper_ci - n_individuals) / n_individuals) * 100
        )

    aic_options <- AICOptions$new(
        model = c("poisson"),
        capture_columns = captures,
    )

    step_options <- StepwiseOptions$new(
        model = c("poisson"),
        capture_columns = captures,
        threshold = c(0.005),
        direction = "both",
        frequency_col_name = "N_ID",
        interaction_limit = 2
    )

    aic_results <- crc(no_cov_sim_data, aic_options)
    step_results <- crc(no_cov_sim_data, step_options)

    selection_df <- data.frame(
        method = c("AIC-Poisson", "Stepwise-Poisson"),
        estimate = c(
            ((aic_results[1, ]$estimate - unknown_pop_size) /
                unknown_pop_size) *
                100,
            ((step_results$estimate - unknown_pop_size) / unknown_pop_size) *
                100
        ),
        lower_ci = c(
            ((aic_results[1, ]$lower_ci - unknown_pop_size) /
                unknown_pop_size) *
                100,
            ((step_results$lower_ci - unknown_pop_size) / unknown_pop_size) *
                100
        ),
        upper_ci = c(
            ((aic_results[1, ]$upper_ci - unknown_pop_size) /
                unknown_pop_size) *
                100,
            ((step_results$upper_ci - unknown_pop_size) / unknown_pop_size) *
                100
        )
    )

    aic_options <- AICOptions$new(
        model = c("negbin"),
        capture_columns = captures,
    )

    step_options <- StepwiseOptions$new(
        model = c("negbin"),
        capture_columns = captures,
        threshold = c(0.005),
        direction = "both",
        frequency_col_name = "N_ID",
        interaction_limit = 2
    )

    aic_nb_results <- crc(no_cov_sim_data, aic_options)
    step_nb_results <- crc(no_cov_sim_data, step_options)

    selection_nb_df <- data.frame(
        method = c("AIC-Negbin", "Stepwise-Negbin"),
        estimate = c(
            ((aic_nb_results[1, ]$estimate - unknown_pop_size) /
                unknown_pop_size) *
                100,
            ((step_nb_results$estimate - unknown_pop_size) / unknown_pop_size) *
                100
        ),
        lower_ci = c(
            ((aic_nb_results[1, ]$lower_ci - unknown_pop_size) /
                unknown_pop_size) *
                100,
            ((step_nb_results$lower_ci - unknown_pop_size) / unknown_pop_size) *
                100
        ),
        upper_ci = c(
            ((aic_nb_results[1, ]$upper_ci - unknown_pop_size) /
                unknown_pop_size) *
                100,
            ((step_nb_results$upper_ci - unknown_pop_size) / unknown_pop_size) *
                100
        )
    )

    combo <- bind_rows(selection_df, selection_nb_df)
    return(bind_rows(drpop_df, combo))
}


set.seed(2024)
n_individuals <- 3e5
n_captures <- 6
p_captures <- c(0.45, 0.35, 0.25, 0.20, 0.4, 0.1)
covariate_df <- data.frame(
    distribution = c("uniform", "uniform"),
    p1 = c(0, 0),
    p2 = c(1, 1),
    dtype = c("integer", "integer")
)

captures <- c(
    "capture_1",
    "capture_2",
    "capture_3",
    "capture_4",
    "capture_5",
    "capture_6"
)

res <- map(
    seq_len(1),
    ~ run_crc_simulations(
        n_individuals = n_individuals,
        n_captures = n_captures,
        p_captures = p_captures,
        covariate_df = covariate_df,
        captures = captures
    )
)

final_results <- bind_rows(res) |>
    group_by(method) |>
    summarise(
        estimate = mean(estimate),
        lower_ci = mean(lower_ci),
        upper_ci = mean(upper_ci)
    )

write.csv(final_results, "final_results.csv", row.names = FALSE)
