################################################################################
# File: simulation_data.R                                                      #
# Project: crcsim                                                              #
# Created Date: 2026-10-05                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-09                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Construct a structured capture-recapture simulation
#'
#' Generate a complete simulated population and return both the observed
#' capture histories and the known all-zero histories used as simulation truth.
#' The data-generating mechanism is homogeneous conditional independence: each
#' capture list has a fixed marginal capture probability and lists are sampled
#' independently. When \\code{alpha} and \\code{beta} are supplied, capture
#' probabilities follow list-specific logistic models conditional on generated
#' covariates.
#'
#' @param n_individuals int: number of individuals to simulate.
#' @param n_captures int: number of binary capture columns (datasets) to
#'   simulate. This must be supplied explicitly for every simulation mode.
#' @param base S3 base capture specification created by
#'   `homogeneous_capture_spec()` or `covariate_capture_spec()`.
#' @param heterogeneity Optional S3 heterogeneity specification created by
#'   `continuous_heterogeneity_spec()` or `latent_class_spec()`.
#' @param dependence Optional S3 dependence specification created by
#'   `pairwise_dependence_spec()` or `graph_dependence_spec()`.
#' @param seed numeric scalar: optional random seed for reproducible output.
#' @details
#' The returned object separates the complete generated population from the
#' observed union and the all-zero population component. The all-zero component
#' is retained as simulation truth and is not available from an observed-only
#' data frame. Fixed-probability simulations record their normalized list
#' probabilities in \code{truth$p_captures}; covariate-dependent simulations
#' record one true probability per individual and capture list in
#' \code{truth$capture_probabilities}. Continuous heterogeneity adds a shared
#' unobserved propensity, while pairwise dependence changes the second list's
#' conditional log odds based on the first list's capture.
#' Latent-class mode assigns each individual to an unobserved class and samples
#' lists independently conditional on that class.
#' @returns A list with `full_data`, `observed_data`, `unobserved_data`,
#'   `observed_frequency_table`, `truth`, and `metadata` components.
#'
#' @export
simulate_data <- function(
    n_individuals,
    n_captures,
    base = homogeneous_capture_spec(),
    heterogeneity = NULL,
    dependence = NULL,
    seed = NULL
) {
    # Validation and Cleaning of Parameters
    validate_and_set_seed(seed)
    validate_n_captures(n_captures)

    if (!inherits(base, "capture_spec")) {
        stop("base must be a capture specification")
    }
    if (inherits(base, "homogeneous_capture_spec")) {
        if (
            !is.null(base$p_captures) && length(base$p_captures) != n_captures
        ) {
            stop("p_captures must have one value per capture list")
        }
        p_captures <- base$p_captures
    } else if (inherits(base, "covariate_capture_spec")) {
        if (length(base$alpha) != n_captures) {
            stop("alpha must have one value per capture list")
        }
        p_captures <- NULL
    } else {
        stop("base must be homogeneous or covariate-dependent")
    }

    latent_class_probabilities <- NULL
    latent_class_p_captures <- NULL
    latent_sd <- 0
    if (!is.null(heterogeneity)) {
        if (inherits(heterogeneity, "continuous_heterogeneity_spec")) {
            latent_sd <- heterogeneity$sd
        } else if (inherits(heterogeneity, "latent_class_spec")) {
            if (ncol(heterogeneity$capture_probabilities) != n_captures) {
                stop(
                    "latent class probabilities need one column per capture list"
                )
            }
            latent_class_probabilities <- heterogeneity$probabilities
            latent_class_p_captures <- heterogeneity$capture_probabilities
        } else {
            stop(
                "heterogeneity must be a supported heterogeneity specification"
            )
        }
    }

    dependence_pair <- NULL
    dependence_strength <- 0
    dependence_edges <- NULL
    if (!is.null(dependence)) {
        if (inherits(dependence, "pairwise_dependence_spec")) {
            dependence_pair <- dependence$pair
            dependence_strength <- dependence$strength
        } else if (inherits(dependence, "graph_dependence_spec")) {
            dependence_edges <- dependence$edges
        } else {
            stop("dependence must be a supported dependence specification")
        }
    }
    validate_dependence_args(
        n_captures = n_captures,
        latent_sd = latent_sd,
        dependence_pair = dependence_pair,
        dependence_strength = dependence_strength,
        dependence_edges = dependence_edges
    )

    # Baseline data (e.g. data that maintains conditional independence)
    if (inherits(base, "covariate_capture_spec")) {
        base_state <- generate_covariate_captures(
            n_individuals = n_individuals,
            n_captures = n_captures,
            p_captures = NULL,
            covariate_ranges = base$covariate_ranges,
            alpha = base$alpha,
            beta = base$beta
        )
    } else {
        base_state <- generate_homogenous_captures(
            n_individuals = n_individuals,
            n_captures = n_captures,
            p_captures = p_captures,
            covariate_ranges = covariate_ranges
        )
    }

    # Extracted default data from the base state
    full_data <- base_state$full_data
    class_assignments <- NULL
    capture_probabilities <- base_state$capture_probabilities
    if (is.null(capture_probabilities)) {
        capture_probabilities <- base_state$base_capture_probabilities
    }
    latent_effects <- rep(0, n_individuals)
    scenario <- base_state$scenario

    capture_columns <- paste0("capture_", seq_len(n_captures))

    # If we want to add unobserved heterogeneity (latent classes)
    latent_class_mode <- !is.null(latent_class_probabilities)
    if (latent_class_mode) {
        validate_latent_class_args(
            latent_sd,
            latent_class_probabilities,
            latent_class_p_captures
        )

        # Normalizing latent class probabilities just in case
        latent_class_probabilities <- latent_class_probabilities /
            sum(latent_class_probabilities)

        latent_state <- apply_latent_classes(
            n_individuals = n_individuals,
            latent_class_probabilities = latent_class_probabilities,
            latent_class_p_captures = latent_class_p_captures
        )
        class_assignments <- latent_state$class_assignments
        base_capture_probabilities <- latent_state$base_capture_probabilities
        latent_state <- apply_latent_effects(
            n_individuals = n_individuals,
            base_capture_probabilities = base_capture_probabilities,
            latent_sd = latent_sd
        )
        capture_probabilities <- latent_state$capture_probabilities
        latent_effects <- latent_state$latent_effects
    } else if (latent_sd > 0) {
        latent_state <- apply_latent_effects(
            n_individuals = n_individuals,
            base_capture_probabilities = capture_probabilities,
            latent_sd = latent_sd
        )
        capture_probabilities <- latent_state$capture_probabilities
        latent_effects <- latent_state$latent_effects
    }

    dependence_edges <- build_dependence_edges(
        dependence_pair = dependence_pair,
        dependence_edges = dependence_edges,
        dependence_strength = dependence_strength
    )

    # If we want to add dependence between captures
    dependence_mode <- !is.null(dependence_pair) || nrow(dependence_edges) > 0
    mechanism_mode <- latent_class_mode || latent_sd > 0 || dependence_mode
    if (mechanism_mode) {
        capture_state <- sample_captures_with_dependence(
            full_data = base_state$full_data,
            n_individuals = n_individuals,
            n_captures = n_captures,
            capture_probabilities = capture_probabilities,
            dependence_edges = dependence_edges
        )
        full_data <- capture_state$full_data
        capture_probabilities <- capture_state$capture_probabilities
    }

    if (mechanism_mode) {
        graph_mode <- is.null(dependence_pair) && nrow(dependence_edges) > 0
        scenario <- get_scenario(
            graph_mode,
            latent_class_mode,
            latent_sd,
            dependence_pair
        )
    }

    observed_data <- extract_captured_data(full_data, capture_columns)
    unobserved_data <- extract_uncaptured_data(full_data, capture_columns)

    result <- assemble_simulation_result(
        full_data = full_data,
        observed_data = observed_data,
        unobserved_data = unobserved_data,
        capture_columns = capture_columns,
        n_captures = n_captures,
        p_captures = base_state$p_captures,
        latent_class_probabilities = latent_class_probabilities,
        latent_class_p_captures = latent_class_p_captures,
        class_assignments = class_assignments,
        alpha = if (inherits(base, "covariate_capture_spec")) {
            base$alpha
        } else {
            NULL
        },
        beta = if (inherits(base, "covariate_capture_spec")) {
            base$beta
        } else {
            NULL
        },
        capture_probabilities = capture_probabilities,
        latent_effects = latent_effects,
        dependence_pair = dependence_pair,
        dependence_strength = dependence_strength,
        dependence_edges = dependence_edges,
        scenario = scenario,
        latent_sd = latent_sd,
        latent_class_mode = latent_class_mode,
        seed = seed
    )
    class(result) <- c("crcsim_simulation", "list")
    return(result)
}

#' Validate the number of capture datasets.
#'
#' @param n_captures Number of capture datasets and capture columns.
#' @return Invisibly returns NULL if the value is valid; otherwise, throws an
#'   error.
#'
#' @keywords internal
#' @noRd
validate_n_captures <- function(n_captures) {
    if (
        is_valid_number(n_captures) ||
            n_captures < 1 ||
            n_captures %% 1 != 0
    ) {
        stop("n_captures must be a positive integer")
    }
    return(invisible(NULL))
}

get_scenario <- function(
    graph_mode,
    latent_class_mode,
    latent_sd,
    dependence_pair
) {
    scenario <- if (graph_mode && latent_class_mode) {
        "latent class heterogeneity with graph dependence"
    } else if (graph_mode) {
        "graph-structured dependence"
    } else if (latent_class_mode) {
        "latent class heterogeneity"
    } else if (latent_sd > 0 && !is.null(dependence_pair)) {
        "latent heterogeneity with direct pairwise dependence"
    } else if (latent_sd > 0) {
        "continuous latent heterogeneity"
    } else {
        "direct pairwise dependence"
    }
    return(scenario)
}

#' Validate and construct dependence edges.
#'
#' @param dependence_pair A pair of dependent captures.
#' @param dependence_edges Existing dependence edges.
#' @param dependence_strength Strength of the dependence.
#' @returns A data frame of validated dependence edges.
#'
#' @keywords internal
#' @noRd
build_dependence_edges <- function(
    dependence_pair,
    dependence_edges,
    dependence_strength
) {
    # Setting up dependence edges if a dependence pair is provided
    if (!is.null(dependence_pair)) {
        dependence_edges <- data.frame(
            from = dependence_pair[[1]],
            to = dependence_pair[[2]],
            strength = dependence_strength
        )
    }
    if (is.null(dependence_edges)) {
        dependence_edges <- data.frame(
            from = integer(),
            to = integer(),
            strength = numeric()
        )
    }
    return(dependence_edges)
}

#' Generate the baseline data. Baseline data assumes that all data is
#' homogenous and maintains conditional independence. It is only accessible
#' when `alpha` and `beta` parameters are not supplied to `simulate_data`.
#'
#' @param n_individuals Number of individuals to generate.
#' @param n_captures Number of capture lists.
#' @param p_captures Fixed list-specific capture probabilities.
#' @param covariate_ranges Covariate-generation specification.
#' @returns A list containing the generated full data, number of captures,
#' fixed capture probabilities, capture probabilities matrix, and scenario
#' label.
#'
#' @keywords internal
#' @noRd
generate_homogenous_captures <- function(
    n_individuals,
    n_captures,
    p_captures,
    covariate_ranges
) {
    if (is.null(n_captures)) {
        stop("n_captures is required without alpha and beta")
    }
    full_data <- create_data(
        n_individuals = n_individuals,
        n_captures = n_captures,
        p_captures = p_captures,
        covariate_ranges = covariate_ranges
    )
    p_captures <- if (is.null(p_captures)) {
        rep(0.5, n_captures)
    } else {
        as.numeric(p_captures)
    }
    return(list(
        full_data = full_data,
        n_captures = n_captures,
        p_captures = p_captures,
        capture_probabilities = NULL,
        base_capture_probabilities = matrix(
            rep(p_captures, each = n_individuals),
            nrow = n_individuals,
            ncol = n_captures
        ),
        scenario = "homogeneous conditional independence"
    ))
}

#' Generate base capture data for a simulation
#'
#' Builds either homogeneous list probabilities or covariate-dependent
#' logistic probabilities before latent heterogeneity or direct dependence is
#' applied.
#'
#' @description This function generates capture data that depends on individual
#' covariates through a logistic model, before any latent heterogeneity or
#' direct dependence is applied. For individual \latex{i}, the probability
#' of capture on list \latex{j} is given by
#' \latex{\text{logit}^{-1}(\alpha_j + \sum_k \beta_{jk} x_{ik})}.
#'
#' @param n_individuals Number of individuals to generate.
#' @param n_captures Number of capture lists supplied by the user.
#' @param p_captures Fixed list-specific capture probabilities.
#' @param covariate_ranges Covariate-generation specification.
#' @param alpha List-specific logistic intercepts.
#' @param beta List-by-covariate logistic coefficient matrix.
#' @param covariate_mode Whether to generate covariate-dependent captures.
#' @returns A list containing the generated data, base probabilities, capture
#'   count, and initial scenario label.
#'
#' @keywords internal
#' @noRd
generate_covariate_captures <- function(
    n_individuals,
    n_captures,
    p_captures,
    covariate_ranges,
    alpha,
    beta
) {
    validate_covariate_capture_model(
        n_captures = n_captures,
        p_captures = p_captures,
        covariate_ranges = covariate_ranges,
        alpha = alpha,
        beta = beta
    )
    covariate_data <- create_data(
        n_individuals = n_individuals,
        n_captures = 1,
        p_captures = 0,
        covariate_ranges = covariate_ranges
    )
    full_data <- covariate_data[
        names(covariate_data) != "capture_1",
        drop = FALSE
    ]
    covariate_columns <- grep(
        "^covariate_",
        names(full_data),
        value = TRUE
    )
    linear_predictors <- sweep(
        as.matrix(full_data[, covariate_columns, drop = FALSE]) %*%
            t(beta),
        2,
        alpha,
        FUN = "+"
    )
    capture_probabilities <- stats::plogis(linear_predictors)
    base_capture_probabilities <- capture_probabilities
    capture_columns <- paste0("capture_", seq_len(n_captures))
    colnames(capture_probabilities) <- capture_columns
    for (capture_column in capture_columns) {
        full_data[[capture_column]] <- as.integer(stats::rbinom(
            n = nrow(full_data),
            size = 1,
            prob = capture_probabilities[, capture_column]
        ))
    }
    full_data <- full_data[
        c(capture_columns, covariate_columns),
        drop = FALSE
    ]
    return(list(
        full_data = full_data,
        n_captures = n_captures,
        p_captures = p_captures,
        capture_probabilities = capture_probabilities,
        base_capture_probabilities = base_capture_probabilities,
        scenario = "covariate-dependent conditional independence"
    ))
}

#' Assign latent classes and their list-specific capture probabilities
#'
#' @param n_individuals Number of individuals to assign to classes.
#' @param latent_class_probabilities Class-assignment probabilities.
#' @param latent_class_p_captures Class-by-list capture-probability matrix.
#' @returns A list containing class assignments and the individual-by-list
#'   base probability matrix.
#'
#' @keywords internal
#' @noRd
apply_latent_classes <- function(
    n_individuals,
    latent_class_probabilities,
    latent_class_p_captures
) {
    class_assignments <- sample(
        seq_along(latent_class_probabilities),
        size = n_individuals,
        replace = TRUE,
        prob = latent_class_probabilities
    )
    return(list(
        class_assignments = class_assignments,
        base_capture_probabilities = latent_class_p_captures[
            class_assignments,
            ,
            drop = FALSE
        ]
    ))
}

#' Apply latent heterogeneity or dependence to base capture probabilities
#'
#' Apply continuous latent heterogeneity to capture probabilities.
#'
#' @param n_individuals Number of individuals to simulate.
#' @param base_capture_probabilities Individual-by-list base probabilities.
#' @param latent_sd Standard deviation of the shared latent propensity.
#' @returns A list containing individual capture probabilities and latent
#'   effects.
#'
#' @keywords internal
#' @noRd
apply_latent_effects <- function(
    n_individuals,
    base_capture_probabilities,
    latent_sd
) {
    if (latent_sd > 0) {
        latent_effects <- stats::rnorm(n_individuals, sd = latent_sd)
        capture_probabilities <- stats::plogis(
            stats::qlogis(base_capture_probabilities) + latent_effects
        )
    } else {
        latent_effects <- rep(0, n_individuals)
        capture_probabilities <- base_capture_probabilities
    }
    return(list(
        capture_probabilities = capture_probabilities,
        latent_effects = latent_effects
    ))
}

#' Sample capture histories using sequential dependence edges.
#'
#' @param full_data Data frame containing non-capture columns and any initial
#'   capture columns.
#' @param n_individuals Number of individuals to simulate.
#' @param n_captures Number of capture lists.
#' @param capture_probabilities Individual-by-list base probabilities.
#' @param dependence_edges Directed dependence graph edges.
#' @returns A list containing the rebuilt data and final capture
#'   probabilities.
#'
#' @keywords internal
#' @noRd
sample_captures_with_dependence <- function(
    full_data,
    n_individuals,
    n_captures,
    capture_probabilities,
    dependence_edges
) {
    capture_columns <- paste0("capture_", seq_len(n_captures))
    full_data <- full_data[
        setdiff(names(full_data), capture_columns),
        drop = FALSE
    ]
    generated_captures <- matrix(
        0L,
        nrow = n_individuals,
        ncol = n_captures
    )
    for (capture_index in seq_len(n_captures)) {
        probabilities <- capture_probabilities[, capture_index]
        incoming_edges <- dependence_edges[
            dependence_edges$to == capture_index,
            ,
            drop = FALSE
        ]
        if (nrow(incoming_edges) > 0) {
            association <- rowSums(vapply(
                seq_len(nrow(incoming_edges)),
                function(edge_index) {
                    return(
                        (generated_captures[,
                            incoming_edges$from[edge_index]
                        ] -
                            0.5) *
                            incoming_edges$strength[edge_index]
                    )
                },
                numeric(n_individuals)
            ))
            probabilities <- stats::plogis(
                stats::qlogis(probabilities) + association
            )
            capture_probabilities[, capture_index] <- probabilities
        }
        generated_captures[, capture_index] <- stats::rbinom(
            n = n_individuals,
            size = 1,
            prob = probabilities
        )
    }
    for (capture_index in seq_len(n_captures)) {
        full_data[[capture_columns[[capture_index]]]] <-
            generated_captures[, capture_index]
    }
    full_data <- full_data[
        c(capture_columns, setdiff(names(full_data), capture_columns)),
        drop = FALSE
    ]
    return(list(
        full_data = full_data,
        capture_probabilities = capture_probabilities
    ))
}

#' Assemble the public simulation result
#'
#' Separates observed and all-zero records and packages data, truth, and
#' metadata into the structure returned by `simulate_data()`.
#'
#' @param full_data Complete generated population.
#' @param observed_data Records captured by at least one list.
#' @param unobserved_data Records with an all-zero capture history.
#' @param capture_columns Names of capture-list columns.
#' @param n_captures Number of capture lists.
#' @param p_captures Fixed list-specific capture probabilities.
#' @param latent_class_probabilities Normalized latent-class probabilities.
#' @param latent_class_p_captures Class-by-list capture-probability matrix.
#' @param class_assignments Individual latent-class assignments.
#' @param alpha List-specific logistic intercepts.
#' @param beta List-by-covariate logistic coefficient matrix.
#' @param capture_probabilities Final individual-by-list probabilities.
#' @param latent_effects Shared latent effects by individual.
#' @param dependence_pair Optional directly dependent pair of lists.
#' @param dependence_strength Direct dependence strength.
#' @param dependence_edges Directed dependence graph edges.
#' @param scenario Resolved scenario label.
#' @param latent_sd Latent propensity standard deviation.
#' @param latent_class_mode Whether latent-class sampling is active.
#' @param seed Random seed used for the simulation.
#' @returns A list containing the simulation views, truth, and metadata.
#'
#' @keywords internal
#' @noRd
assemble_simulation_result <- function(
    full_data,
    observed_data,
    unobserved_data,
    capture_columns,
    n_captures,
    p_captures,
    latent_class_probabilities,
    latent_class_p_captures,
    class_assignments,
    alpha,
    beta,
    capture_probabilities,
    latent_effects,
    dependence_pair,
    dependence_strength,
    dependence_edges,
    scenario,
    latent_sd,
    latent_class_mode,
    seed
) {
    return(list(
        full_data = full_data,
        observed_data = observed_data,
        unobserved_data = unobserved_data,
        observed_frequency_table = build_contingency_table(observed_data),
        truth = list(
            n_total = as.integer(nrow(full_data)),
            n_observed = as.integer(nrow(observed_data)),
            n_unobserved = as.integer(nrow(unobserved_data)),
            p_captures = p_captures,
            latent_class_probabilities = latent_class_probabilities,
            latent_class_p_captures = latent_class_p_captures,
            class_assignments = class_assignments,
            alpha = alpha,
            beta = beta,
            capture_probabilities = capture_probabilities,
            latent_effects = latent_effects,
            dependence_pair = dependence_pair,
            dependence_strength = dependence_strength,
            dependence_edges = dependence_edges
        ),
        metadata = list(
            capture_columns = capture_columns,
            n_captures = as.integer(n_captures),
            scenario = scenario,
            latent_sd = latent_sd,
            latent_class_count = if (latent_class_mode) {
                nrow(latent_class_p_captures)
            } else {
                NULL
            },
            seed = seed
        )
    ))
}

validate_covariate_dependent_args <- function(alpha, beta) {
    if (!is.null(alpha) || !is.null(beta)) {
        stop("latent classes cannot be combined with alpha and beta")
    }
    return(invisible(NULL))
}

#' Validate the dependence arguments
#'
#' @param n_captures Number of capture lists. Checks it is a positive integer.
#' @param latent_sd Standard deviation of the latent effects. Checks it is a
#' non-negative numeric scalar.
#' @param dependence_pair Ordered pair of dependent capture lists.
#' @param dependence_strength Strength of the pairwise dependence.
#' @param dependence_edges Data frame specifying multiple dependence edges.
#' @return Invisibly returns NULL if the arguments are valid; otherwise, throws
#' an error.
#'
#' @noRd
#' @keywords internal
validate_dependence_args <- function(
    n_captures,
    latent_sd,
    dependence_pair,
    dependence_strength,
    dependence_edges
) {
    if (is_valid_number(latent_sd) || latent_sd < 0) {
        stop("latent_sd must be a non-negative numeric scalar")
    }
    if (is_valid_number(dependence_strength)) {
        stop("dependence_strength must be a numeric scalar")
    }
    if (!is.null(dependence_pair) && !is.null(dependence_edges)) {
        stop("supply either dependence_pair or dependence_edges, not both")
    }
    if (is.null(dependence_pair) && is.null(dependence_edges)) {
        return(invisible(NULL))
    }
    if (!is.null(dependence_edges)) {
        validate_dependence_edges(dependence_edges)
        if (!is.null(n_captures) && any(dependence_edges$to > n_captures)) {
            stop("dependence_edges must refer to existing capture lists")
        }
        return(invisible(NULL))
    }
    validate_dependence_pair(dependence_pair)
    if (!is.null(n_captures) && dependence_pair[2] > n_captures) {
        stop("dependence_pair must refer to existing capture lists")
    }
    return(invisible(NULL))
}

#' Validates the dependence edges has the correct structure and values. This
#' means it has more than one row, contains numeric columns for from, to, and
#' strength, and that the from and to columns represent ordered graph edges.
#'
#' @param dependence_edges A data frame containing the dependence edges to
#' validate.
#' @return Invisibly returns NULL if the dependence edges are valid.
#' @throws An error if the dependence edges are not valid.
#'
#' @noRd
#' @keywords internal
validate_dependence_edges <- function(dependence_edges) {
    if (
        !is.data.frame(dependence_edges) ||
            !all(c("from", "to", "strength") %in% names(dependence_edges))
    ) {
        stop("dependence_edges must contain from, to, and strength columns")
    }
    if (
        nrow(dependence_edges) < 1 ||
            !is.numeric(dependence_edges$from) ||
            !is.numeric(dependence_edges$to) ||
            !is.numeric(dependence_edges$strength) ||
            anyNA(dependence_edges) ||
            any(dependence_edges$from %% 1 != 0) ||
            any(dependence_edges$to %% 1 != 0) ||
            any(dependence_edges$from < 1) ||
            any(dependence_edges$to < 1) ||
            any(dependence_edges$from >= dependence_edges$to)
    ) {
        stop("dependence_edges must contain ordered numeric graph edges")
    }
    return(invisible(NULL))
}

#' Validates that the dependence pair has the correct structure and values. This
#' means it is a numeric vector of length 2 with positive integer values and
#' that the first element is less than the second element.
#'
#' @param dependence_pair A numeric vector containing the dependence pair to
#' validate.
#' @return Invisibly returns NULL if the dependence pair is valid.
#' @throws An error if the dependence pair is not valid.
#'
#' @noRd
#' @keywords internal
validate_dependence_pair <- function(dependence_pair) {
    if (
        !is.numeric(dependence_pair) ||
            length(dependence_pair) != 2 ||
            anyNA(dependence_pair) ||
            any(dependence_pair %% 1 != 0) ||
            any(dependence_pair < 1) ||
            dependence_pair[1] >= dependence_pair[2]
    ) {
        stop("dependence_pair must be an ordered pair of distinct list indices")
    }
    return(invisible(NULL))
}

#' Validates the arguments for latent class models. Ensures that the latent
#' class probabilities and capture probabilities are correctly specified and
#' consistent with the number of captures.
#'
#' @param n_captures The number of capture lists.
#' @param p_captures The capture probabilities for each list (ignored if latent
#' classes are used).
#' @param alpha The alpha parameter for the beta distribution (ignored if
#' latent classes are used).
#' @param beta The beta parameter for the beta distribution (ignored if latent
#' classes are used).
#' @param latent_class_probabilities A numeric vector of latent class
#' probabilities.
#' @param latent_class_p_captures A matrix of capture probabilities for each
#' latent class.
#' @return Invisibly returns NULL if the arguments are valid.
#' @throws An error if the arguments are not valid.
#'
#' @noRd
#' @keywords internal
validate_latent_class_args <- function(
    latent_sd,
    latent_class_probabilities,
    latent_class_p_captures
) {
    if (!is_valid_number(latent_sd) || latent_sd < 0) {
        stop("latent_sd must be a non-negative number")
    }

    # Check the sizing and existance of latent class matrices
    if (
        is.null(latent_class_probabilities) && is.null(latent_class_p_captures)
    ) {
        return(invisible(NULL))
    } else if (is.null(latent_class_probabilities)) {
        stop(
            "latent_class_p_captures cannot be supplied without ",
            "latent_class_probabilities"
        )
    } else if (is.null(latent_class_p_captures)) {
        stop(
            "if latent_class_probabilities are supplied, we must supply ",
            "latent_class_p_captures."
        )
    } else if (
        nrow(latent_class_p_captures) != length(latent_class_probabilities)
    ) {
        stop("latent_class_p_captures needs one row per latent class")
    }

    # Validate values of latent class matrices
    if (
        !is.numeric(latent_class_probabilities) ||
            length(latent_class_probabilities) < 2 ||
            anyNA(latent_class_probabilities) ||
            any(latent_class_probabilities < 0) ||
            sum(latent_class_probabilities) <= 0
    ) {
        stop("latent_class_probabilities must be non-negative probabilities")
    }
    if (
        !is.numeric(latent_class_p_captures) ||
            !is.matrix(latent_class_p_captures) ||
            anyNA(latent_class_p_captures) ||
            any(latent_class_p_captures < 0 | latent_class_p_captures > 1)
    ) {
        stop("latent_class_p_captures must be a probability matrix")
    }
    return(invisible(NULL))
}

#' Validate the seed provided
#'
#' @param seed An integer scalar or NULL.
#' @return Invisibly returns NULL if the seed is valid; otherwise, throws an
#' error.
#'
#' @keywords internal
#' @noRd
validate_and_set_seed <- function(seed) {
    if (
        !is.null(seed) &&
            (!is.numeric(seed) ||
                length(seed) != 1 ||
                is.na(seed) ||
                seed %% 1 != 0)
    ) {
        stop("seed must be NULL or an integer scalar")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    return(invisible(NULL))
}

validate_covariate_capture_model <- function(
    n_captures,
    p_captures,
    covariate_ranges,
    alpha,
    beta
) {
    if (is.null(alpha) || is.null(beta)) {
        stop("alpha and beta must be supplied together")
    }
    if (!is.null(p_captures)) {
        stop("p_captures cannot be supplied with alpha and beta")
    }
    if (is.null(covariate_ranges)) {
        stop("covariate_ranges is required with alpha and beta")
    }
    if (!is.numeric(alpha) || length(alpha) == 0 || anyNA(alpha)) {
        stop("alpha must be a non-empty numeric vector without NA")
    }
    if (!is.numeric(beta) || !is.matrix(beta) || anyNA(beta)) {
        stop("beta must be a numeric matrix without NA")
    }
    if (nrow(beta) != length(alpha)) {
        stop("beta must have one row for each alpha value")
    }
    if (ncol(beta) != nrow(covariate_ranges)) {
        stop("beta must have one column for each covariate")
    }
    if (is_valid_number(n_captures) || n_captures != length(alpha)) {
        stop("n_captures must be equal to the length of alpha")
    }
    return(invisible(NULL))
}
