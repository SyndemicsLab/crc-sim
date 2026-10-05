################################################################################
# File: simulate.R                                                             #
# Project: crcsim                                                              #
# Created Date: 2026-08-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-09-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Run the CRC Simulations for the provided number of individuals, covariates,
#' and captures.
#'
#' @param n_individuals The number of individuals to simulate.
#' @param n_covariates The covariates that accompany the captures.
#' @param captures the captures that accompany each covariate.
#' @return The results of the CRC Simulation
#'
#' @importFrom dplyr select rename mutate starts_with
#' @export
simulate <- function(
    n_individuals,
    n_covariates,
    captures,
    link_function = "logit",
    cov_func = covariate_conditions
) {
    cov_cols <- paste0("covariate_", seq_len(n_covariates))
    sim_data <- build_covariate_data(n_individuals, n_covariates) |>
        captures_from_covariates(length(captures), cov_func) |>
        relocate_covariates(cov_cols) |>
        all_int_cols_to_numeric() |>
        extract_captured_data(captures)

    no_cov_sim_data <- sim_data |>
        select(-starts_with("covariate_")) |>
        build_contingency_table() |>
        suppress_data(10)

    known_pop_size <- nrow(sim_data)
    unknown_pop_size <- n_individuals - known_pop_size

    internal_estimates <- drpop::popsize(
        sim_data,
        K = length(captures),
        funcname = "logit",
        nfolds = 2,
        margin = 0.02
    )$result |>
        group_by(method) |>
        summarise(
            estimate = mean(n),
            lower_ci = min(cin.l),
            upper_ci = max(cin.u)
        )

    # internal_estimates <- lapply(
    #     c("plugin", "doubly_robust", "tmle"),
    #     function(method) {
    #         estimate <- estimate_capture_prob(
    #             data = sim_data,
    #             n_lists = length(captures),
    #             method = method,
    #             func = link_function,
    #             nfolds = 2,
    #             margin = 0.005,
    #             seed = 1
    #         )
    #         return(estimate)
    #     }
    # )
    # internal_estimates <- bind_rows(internal_estimates) |>
    #     mutate(
    #         method = rep(
    #             c("PI", "DR", "TMLE"),
    #             vapply(internal_estimates, nrow, integer(1))
    #         ),
    #         lower_ci = ci_l,
    #         upper_ci = ci_u
    #     ) |>
    #     select(method, n, lower_ci, upper_ci) |>
    #     rename(
    #         estimate = n
    #     )

    internal_df <- internal_estimates |>
        mutate(
            estimate = ((estimate - n_individuals) / n_individuals) * 100,
            lower_ci = ((lower_ci - n_individuals) / n_individuals) * 100,
            upper_ci = ((upper_ci - n_individuals) / n_individuals) * 100
        )

    aic_options <- LoglinearOptions$new(
        capture_columns = captures,
        model_family = "poisson",
        selection_method = "aic"
    )

    step_options <- LoglinearOptions$new(
        capture_columns = captures,
        model_family = "poisson",
        selection_method = "stepwise",
        selection_options = list(
            direction = "both",
            interaction_limit = 2
        )
    )

    aic_results <- crc(no_cov_sim_data, aic_options)
    step_results <- crc(no_cov_sim_data, step_options)

    # nolint start: indentation_linter
    # Turning off indentation lint because the double parentheses are causing
    # issues. The code is correct and functions as expected.
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
    # nolint end: indentation_linter

    aic_options <- LoglinearOptions$new(
        capture_columns = captures,
        model_family = "negbin",
        selection_method = "aic"
    )

    step_options <- LoglinearOptions$new(
        capture_columns = captures,
        model_family = "negbin",
        selection_method = "stepwise",
        selection_options = list(
            direction = "both",
            interaction_limit = 2
        )
    )

    aic_nb_results <- crc(no_cov_sim_data, aic_options)
    step_nb_results <- crc(no_cov_sim_data, step_options)

    # nolint start: indentation_linter
    # Turning off indentation lint because the double parentheses are causing
    # issues. The code is correct and functions as expected.
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
    # nolint end: indentation_linter

    combo <- bind_rows(selection_df, selection_nb_df)
    return(bind_rows(internal_df, combo))
}
