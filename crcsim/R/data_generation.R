################################################################################
# File: data_generation.R                                                      #
# Project: crcsim                                                              #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-22                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' @title Create the data for a Scenario
#'
#' @description This function is the core entry for generating
#' capture-recapture data based on a specified simulation scenario, including
#' base, latent, and dependence components.
#'
#' @details The function allows the creation of multiple different data
#' scenarios. Each scenario is composed of a "base" component and then
#' additional dependency components. The "base" component represents the
#' fundamental capture probabilities, which can be homogeneous or
#' covariate-dependent. The dependency components are represented through
#' heterogeneityor a direct graph structure. The resulting data includes
#' the generated capture histories and the named scenario. All scenarios are
#' governed through the `scenario` S3 class.
#'
#' @param n_pop integer: The number of individuals in the population.
#' @param n_capture_cols integer: The number of capture lists.
#' @param scenario A simulation scenario object created by
#' `simulation_scenario()`.
#' @param seed Optional random seed for reproducibility.
#' @returns A list containing the generated capture histories and any
#' additional scenario-specific data.
#'
#' @export
create_scenario_data <- function(
    n_pop,
    n_capture_cols,
    scenario = simulation_scenario(),
    seed = NULL
) {
    # Validation and Cleaning of General Parameters
    validate_and_set_seed(seed)
    validate_n_pop(n_pop)
    validate_n_capture_cols(n_capture_cols)
    data_list <- list()
    base_scenario <- scenario$base_scenario
    if (inherits(base_scenario, "homogeneous_scenario")) {
        data_list <- create_homogeneous_scenario(
            n_pop,
            n_capture_cols,
            base_scenario$capture_probs,
            n_covariates = base_scenario$n_covariates
        )
    } else if (inherits(base_scenario, "covariate_scenario")) {
        data_list <- create_covariate_scenario(
            n_pop,
            n_capture_cols,
            base_scenario$cov_dists,
            base_scenario$alpha,
            base_scenario$beta
        )
    } else {
        stop("Unsupported base scenario")
    }

    if (!is.null(scenario$heterogeneity_scenario)) {
        # Do Heterogeneity Details
        data_list <- add_latent_scenario(
            data_list,
            scenario$heterogeneity_scenario
        )
    }

    if (!is.null(scenario$dependence_scenario)) {
        # Do Dependence Details
        data_list <- add_graph_dependence_scenario(
            data_list,
            scenario$dependence_scenario
        )
    }
    return(data_list)
}

#' @title Create a Homogeneous Capture Scenario
#'
#' @description This function creates a simulated dataset under a homogeneous
#' capture scenario, where each capture column has the same probability of
#' capture for all individuals.
#'
#' @details The homogeneous capture scenario assumes that the probability of
#' capture is the same for all individuals across each capture column.
#' Covariates, if specified, do not influence the capture probabilities.
#'
#' @param n_pop integer: The number of people to simulate
#' @param n_capture_cols integer: The number of binary capture columns to
#' simulate
#' @param capture_probs numeric vector: The capture probabilities for each
#' capture column
#' @param cov_dists data.frame: Optional covariate specification table
#' @param n_covariates integer: The number of covariates to generate if
#' \code{cov_dists} is not provided
#' @return list: A list containing the simulated dataset and the scenario type
#'
#' @keywords internal
create_homogeneous_scenario <- function(
    n_pop,
    n_capture_cols,
    capture_probs,
    cov_dists = NULL,
    n_covariates = 0
) {
    if (length(capture_probs) == 1) {
        capture_probs <- rep(capture_probs, n_capture_cols)
    }
    if (n_covariates > 0 && is.null(cov_dists)) {
        cov_dists <- data.frame(
            distribution = rep("uniform", n_covariates),
            p1 = rep(0.0, n_covariates),
            p2 = rep(1.0, n_covariates),
            dtype = rep("integer", n_covariates)
        )
    }

    full_data <- create_dataset(n_pop, n_capture_cols, capture_probs, cov_dists)
    return(list(
        full_data = full_data,
        scenario = "homogeneous"
    ))
}

#' @title Create a Covariate Scenario
#'
#' @description This function generates capture data that depends on individual
#' covariates through a logistic model, before any latent heterogeneity or
#' direct dependence is applied. For individual \latex{i}, the probability
#' of capture on list \latex{j} is given by
#' \latex{\text{logit}^{-1}(\alpha_j + \sum_k \beta_{jk} x_{ik})}.
#'
#' @param n_pop Number of individuals to generate.
#' @param n_captures_cols Number of capture lists supplied by the user.
#' @param cov_dists Covariate-generation specification.
#' @param alpha List-specific logistic intercepts.
#' @param beta List-by-covariate logistic coefficient matrix.
#' @returns A list containing the generated data and initial scenario label.
#'
#' @keywords internal
create_covariate_scenario <- function(
    n_pop,
    n_capture_cols,
    cov_dists,
    alpha,
    beta
) {
    validate_covariate_scenario(n_capture_cols, cov_dists, alpha, beta)
    capture_columns <- paste0("capture_", seq_len(n_captures))

    cov_data <- create_dataset(n_pop, 1, 0, cov_dists) |>
        select(-capture_1)
    full_data <- cov_data

    covariate_columns <- colnames(cov_data)

    # This calculates the covariate data * beta coefficients and then applies
    # the alpha intercept across each column (the 2 dim) of the matrix.
    # sweep is just a poorly named function to apply a vector across a matrix.
    capture_probabilities <- sweep(
        as.matrix(cov_data) %*% t(beta),
        2,
        alpha,
        FUN = "+"
    ) |>
        stats::plogis()

    colnames(capture_probabilities) <- capture_columns

    full_data[capture_columns] <- matrix(
        as.integer(stats::rbinom(
            n = length(capture_probabilities),
            size = 1,
            prob = capture_probabilities
        )),
        nrow = nrow(capture_probabilities),
        dimnames = dimnames(capture_probabilities)
    ) |>
        as.data.frame(capture_probabilities)

    full_data <- dplyr::relocate(
        full_data,
        any_of(covariate_columns),
        .after = last_col()
    )

    return(list(
        full_data = full_data,
        scenario = "covariate-dependent conditional independence"
    ))
}

#' @title Create Dataset
#'
#' @description Data generation tool. This creates a simulated dataset for
#' testing various Capture-Recapture models.
#'
#' @details The user is able to specify the number of individuals, the number
#' of captures, and optional numerical covariates generated from per-covariate
#' distribution settings. The returned value represents the "raw" simulated
#' data where each row is a single sample. Frequency tables can also be
#' generated by setting \code{return_frequency_table} to \code{TRUE}, which
#' will aggregate identical rows and append a frequency column \code{N_ID}.
#'
#' @param n_pop int: number of participants to simulate
#' @param n_capture_cols int: number of binary capture columns to simulate
#' @param capture_probs numeric vector: optional capture probabilities of length
#'   \code{n_capture_cols}. Defaults to 0.5 for each capture column.
#' @param cov_dists data.frame: optional covariate specification table
#'   with columns \code{distribution}, \code{p1}, \code{p2}, and \code{dtype}.
#'   Each row defines one numerical covariate named \code{covariate_i}.
#'   Supported distributions are \code{uniform} and \code{normal}. For
#'   \code{uniform}, \code{p1 = low} and \code{p2 = high}. For \code{normal},
#'   \code{p1 = mean} and \code{p2 = standard deviation}.
#' @param return_frequency_table logical: when \code{TRUE}, aggregates identical
#'   rows and appends frequency column \code{N_ID}. Default is \code{FALSE}.
#'
#' @returns a tibble
#'
#' @export
create_dataset <- function(
    n_pop,
    n_capture_cols,
    capture_probs = NULL,
    cov_dists = NULL,
    freq_table = FALSE
) {
    validate_create_data_params(n_pop, n_capture_cols, cov_dists, freq_table)
    capture_probs <- validate_capture_probs(capture_probs, n_capture_cols)
    cov_dists <- validate_cov_dists(cov_dists)

    records <- build_records(n_pop, capture_probs, cov_dists)

    if (isTRUE(freq_table)) {
        return(create_contingency_table(records))
    }

    return(records)
}

#' @title Create an Alpha Intercept Vector
#'
#' @description Create an intercept vector to be used with covariant
#' relationships to capture probabilities.
#'
#' @details The alpha vector returned has length `n_capture_cols` and is filled
#' with 0.0. This is because we assume the intercept for each capture column is
#' initially zero. Users are intended to use this as a shell vector and modify
#' to their scenario.
#'
#' @param n_capture_cols Number of capture columns.
#' @param random Whether to generate random intercepts instead of zeros.
#' Default is \code{FALSE}.
#' @returns A numeric vector of zeros with length `n_capture_cols`.
#'
#' @export
create_alpha <- function(n_capture_cols, random = FALSE) {
    if (isTRUE(random)) {
        return(stats::runif(n_capture_cols))
    }
    return(rep(0.0, n_capture_cols))
}

#' @title Create a Beta Matrix
#'
#' @description Create a beta coefficient matrix to establish covariant
#' relationships with capture probabilities.
#'
#' @details The beta matrix returned has dimensions `n_capture_cols` by
#' `n_covariates` and is filled with 0.0. This is because we assume each
#' covariate has no effect on the capture probabilities. Users are intended to
#' use this as a shell matrix and modify to their scenario.
#'
#' @param n_capture_cols Number of capture columns.
#' @param n_covariates Number of covariates.
#' @param random Whether to generate random coefficients instead of zeros.
#' Default is \code{FALSE}.
#' @returns A matrix of zeros with dimensions `n_capture_cols` by
#' `n_covariates`.
#'
#' @export
create_beta <- function(n_capture_cols, n_covariates, random = FALSE) {
    if (isTRUE(random)) {
        return(matrix(
            stats::runif(n_capture_cols * n_covariates),
            nrow = n_capture_cols,
            ncol = n_covariates
        ))
    }
    return(matrix(0.0, nrow = n_capture_cols, ncol = n_covariates))
}

#' @title Build Records.
#'
#' @description This builds the full set of records of `n_pop` people and their
#' capture history and the corresponding covariates.
#'
#' @details This is a wrapper function to sample and combine `n_pop` individual
#' records. Each record is a row in the final dataframe and consists of
#' \code{ncols(capture_probs) + \nrows(cov_dists)} columns. This is because a
#' record is identified by each unique capture column and each covariate column.
#' This function is traditionally wrapped by \code{create_data} but exposed for
#' users who want more granular control over building the simulation data.
#' Importantly, the capture_probs and covariate_specs are not validated so
#' errors may not be as descriptive.
#'
#' @param n_pop int: number of participants to simulate
#' @param capture_probs numeric vector: capture probabilities of length equal to
#'  the number of capture columns
#' @param covariate_specs data.frame: covariate specification table with one row
#' per covariate and columns \code{distribution}, \code{p1}, \code{p2}, and
#' \code{dtype}.
#' @returns a tibble where each row is an individual record
#'
#' @keywords internal
#' @export
build_records <- function(n_pop, capture_probs, cov_dists) {
    records <- purrr::map(
        seq_len(n_pop),
        create_single_record,
        capture_probs = capture_probs,
        cov_dists = cov_dists
    ) |>
        purrr::list_rbind()

    return(records)
}

#' Create One Individual Simulation Record. This simulates a single record of
#' data for one individual, including the capture indicators and numerical
#' covariates. The return type is a tibble row.
#' This is an internal function used within \code{build_records}.
#'
#' @keywords internal
#' @noRd
create_single_record <- function(index, capture_probs, cov_dists) {
    record_values <- capture_probs |>
        binom_sample(prefix = "capture_") |>
        as.list()

    covariate_values <- dist_df_sample(cov_dists, prefix = "covariate_")

    return(invisible(tibble::as_tibble_row(c(record_values, covariate_values))))
}

#' @title Create Contingency Table
#'
#' @description A function to aggregate individual records into a contingency
#' table. This function takes the row level records produced by
#' \code{build_records} and aggregates identical rows while appending a
#' frequency column \code{N_ID}.
#'
#' @param records tibble: raw data where each row is an individual record with
#' capture and covariate category columns.
#' @return a contingency table tibble where identical rows in the input are
#' aggregated and a frequency column \code{N_ID} is appended.
#'
#' @importFrom dplyr count across everything
#'
#' @keywords internal
#' @export
create_contingency_table <- function(records) {
    return(count(records, across(everything()), name = "N_ID"))
}


#' Sample Captures Based on Covariates. This function generates a data frame of
#' simulated capture data for a specified number of individuals based on the
#' provided covariate data. Each row represents an individual, and each column
#' represents a capture event. The values in the data frame are binary (0 or 1),
#' indicating whether an individual was captured (1) or not (0) for each
#' capture event. The capture probabilities are determined by the covariate
#' values, specifically the "sex" column, where the probability of capture is
#' higher for one sex compared to the other.
#'
#' @param data_table data.frame: A data frame containing covariate data for each
#' individual. The data frame must contain a "sex" column with binary values (0
#' or 1) indicating the sex of each individual.
#' @param n_cols int: The number of capture columns to generate. Each column
#' represents a capture event.
#' @param covariate_conditions function: function that takes \code{data_table}
#' and returns a numeric probability vector \code{capture_prob} in [0, 1].
#' @return A data frame where each row represents an individual and each column
#' represents a capture event. The values are binary (0 or 1), indicating
#' whether an individual was captured (1) or not (0) for each capture event,
#' based on the covariate values in the input data frame.
#'
#' @importFrom stats rbinom
#' @export
captures_from_covariates <- function(
    data_table,
    n_cols,
    cov_func
) {
    if (!is.data.frame(data_table)) {
        stop("data_table must be a data.frame")
    }
    if (!is.numeric(n_cols) || length(n_cols) != 1 || n_cols %% 1 != 0) {
        stop("n_cols must be an integer scalar")
    }
    if (!is.function(cov_func)) {
        stop("cov_func must be a function")
    }
    if (n_cols <= 0) {
        return(data_table)
    }

    capture_prob <- cov_func(data_table)
    if (!is.numeric(capture_prob) || length(capture_prob) != nrow(data_table)) {
        stop(paste0(
            "cov_func must return a numeric vector",
            " of length nrow(data_table)"
        ))
    }
    if (any(is.na(capture_prob) | capture_prob < 0 | capture_prob > 1)) {
        stop(paste0(
            "cov_func must return probabilities between",
            " 0 and 1 without NA"
        ))
    }

    capture_cols <- paste0("capture_", seq_len(n_cols))
    for (col in capture_cols) {
        if (!col %in% colnames(data_table)) {
            data_table[[col]] <- as.integer(rbinom(
                n = nrow(data_table),
                size = 1,
                prob = capture_prob
            ))
        }
    }
    return(data_table)
}

#' Default Covariate Conditions for Capture Probability.
#' This function computes \code{capture_prob} from four covariates using a
#' logistic transform.
#'
#' @param data_table data.frame containing columns \code{covariate_1},
#' \code{covariate_2}, \code{covariate_3}, and \code{covariate_4}.
#' @param intercept The intercept value to update the covariates.
#' @return numeric vector \code{capture_prob} with values in (0, 1).
#' @export
covariate_conditions <- function(data_table, intercept = 1) {
    z1_values <- data_table[["covariate_1"]]
    z2_values <- data_table[["covariate_2"]]
    z3_values <- data_table[["covariate_3"]]
    z4_values <- data_table[["covariate_4"]]

    v <- exp(
        intercept +
            0.5 * z1_values +
            0.5 * z2_values -
            0.5 * z3_values -
            0.5 * z4_values
    )

    capture_prob <- v / (1 + v)
    return(capture_prob)
}

#' A function used to check if a column in a data frame is binary
#' (i.e., contains only 0s and 1s).
#'
#' @param column a vector representing a column in a data frame
#' @returns a boolean indicating whether the column is binary
#'
#' @keywords internal
#' @export
is_column_binary <- function(column) {
    unique_values <- unique(column)
    return(length(unique_values) == 2 && all(unique_values %in% c(0, 1)))
}

#' A function used to check if the provided data is a frequency table. This is
#' necessary as the CRC estimation functions require a frequency table as input.
#'
#' @param data a data frame
#' @param frequency_column a string specifying the name of the frequency column
#' in the data frame
#' @returns a boolean indicating whether the data frame is a frequency table
#' (i.e., contains the specified frequency column and that column is numeric)
#'
#' @keywords internal
#' @export
is_frequency_table <- function(data, frequency_column) {
    if (
        frequency_column %in%
            names(data) &&
            is.numeric(data[[frequency_column]])
    ) {
        return(TRUE)
    }
    return(FALSE)
}
