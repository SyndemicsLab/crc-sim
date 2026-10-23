################################################################################
# File: simulation_data.R                                                      #
# Project: crcsim                                                              #
# Created Date: 2026-10-05                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-22                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Generate a full capture-recapture simulation data frame
#'
#' Generate one row per individual with capture histories and generated
#' covariates. The data-generating mechanism is homogeneous conditional
#' independence by default: each capture list has a fixed marginal capture
#' probability and lists are sampled independently. When list-specific
#' covariate coefficients are supplied, capture probabilities follow logistic
#' models conditional on generated covariates.
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
#' The returned data frame is the complete generated population and includes
#' all-zero capture histories. Observed and unobserved views should be derived
#' by downstream analysis code.
#' @returns A data frame with one row per simulated individual.
#'
#' @noRd
simulate_data <- function(
    n_individuals,
    n_captures,
    n_covariates = 0,
    base = homogeneous_capture_spec(),
    heterogeneity = NULL,
    dependence = NULL,
    seed = NULL
) {
    # Validation and Cleaning of General Parameters
    validate_and_set_seed(seed)
    validate_n_individuals(n_individuals)
    validate_n_captures(n_captures)

    p_captures <- validate_base(base, n_captures)
    dependence_details <- validate_dependence(dependence)

    # Baseline data
    if (inherits(base, "covariate_capture_spec")) {
        if (n_covariates != nrow(base$covariate_ranges)) {
            stop(paste(
                "The number of covariates specified does not match the",
                "provided covariate ranges!!"
            ))
        }
        base_state <- generate_covariate_captures(
            n_pop = n_individuals,
            n_captures = n_captures,
            p_captures = NULL,
            cov_dists = base$covariate_ranges,
            alpha = base$alpha,
            beta = base$beta
        )
    } else {
        base_state <- create_dataset(
            n_individuals = n_individuals,
            n_captures = n_captures,
            p_captures = p_captures,
            covariate_ranges = data.frame(
                distribution = rep("uniform", n_covariates),
                p1 = rep(0.0, n_covariates),
                p2 = rep(1.0, n_covariates),
                dtype = rep("integer", n_covariates)
            )
        )
    }

    # Extract the generated data and individual capture probabilities.
    full_data <- base_state$full_data
    capture_probabilities <- base_state$capture_probabilities %||%
        base_state$base_capture_probabilities

    if (!is.null(heterogeneity)) {
        latent_details <- build_latent_details(
            heterogeneity,
            n_captures,
            n_individuals,
            capture_probabilities
        )
        capture_probabilities <- latent_details$capture_probabilities %||%
            capture_probabilities
    }

    dependence_edges <- build_dependence_edges(
        dependence_pair = dependence_details$dependence_pair,
        dependence_edges = dependence_details$dependence_edges,
        dependence_strength = dependence_details$dependence_strength
    )

    if (!is.null(heterogeneity) || !is.null(dependence)) {
        capture_state <- sample_captures_with_dependence(
            full_data = full_data,
            n_individuals = n_individuals,
            n_captures = n_captures,
            capture_probabilities = capture_probabilities,
            dependence_edges = dependence_edges
        )
        full_data <- capture_state$full_data
    }

    return(full_data)
}

#' Validate the heterogeneity specification.
#'
#' @param heterogeneity A heterogeneity specification object.
#' @param n_captures Number of capture datasets and capture columns.
#' @return A list containing validated latent class probabilities, latent class
#' capture probabilities, and latent standard deviation.
#'
#' @keywords internal
#' @noRd
validate_heterogeneity <- function(heterogeneity, n_captures) {
    latent_class_probabilities <- NULL
    latent_class_p_captures <- NULL
    latent_sd <- 0
    if (!is.null(heterogeneity)) {
        if (inherits(heterogeneity, "continuous_heterogeneity_spec")) {
            latent_sd <- heterogeneity$sd
        } else if (inherits(heterogeneity, "latent_class_spec")) {
            if (ncol(heterogeneity$capture_probabilities) != n_captures) {
                stop(
                    paste(
                        "latent class probabilities need one column per",
                        "capture list"
                    )
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
    return(list(
        latent_class_probabilities = latent_class_probabilities,
        latent_class_p_captures = latent_class_p_captures,
        latent_sd = latent_sd
    ))
}

#' Validate the dependence specification.
#'
#' @param dependence A dependence specification object.
#' @return A list containing validated dependence pair, dependence strength,
#' and dependence edges.
#'
#' @keywords internal
#' @noRd
validate_dependence <- function(dependence) {
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
    return(list(
        dependence_pair = dependence_pair,
        dependence_strength = dependence_strength,
        dependence_edges = dependence_edges
    ))
}

#' Build latent details for the simulation.
#'
#' @param heterogeneity A heterogeneity specification object.
#' @param n_captures Number of capture datasets and capture columns.
#' @param n_individuals Number of individuals in the simulation.
#' @param capture_probabilities Base capture probabilities for each individual.
#' @return A list containing class assignments, latent effects, and updated
#' capture probabilities.
#'
#' @keywords internal
#' @noRd
build_latent_details <- function(
    heterogeneity,
    n_captures,
    n_individuals,
    capture_probabilities
) {
    class_assignments <- NULL
    latent_effects <- NULL

    het_details <- validate_heterogeneity(heterogeneity, n_captures)

    # If we want to add unobserved heterogeneity (latent classes)
    if (!is.null(het_details$latent_class_probabilities)) {
        lsd <- het_details$latent_sd
        lcp <- het_details$latent_class_probabilities
        lpc <- het_details$latent_class_p_captures
        validate_latent_class_args(lsd, lcp, lpc)

        # Normalizing latent class probabilities just in case
        lcp <- lcp / sum(lcp)

        l_classes <- apply_latent_classes(
            n_individuals = n_individuals,
            latent_class_probabilities = lcp,
            latent_class_p_captures = lpc
        )
        class_assignments <- l_classes$class_assignments
        base_capture_probabilities <- l_classes$base_capture_probabilities
        l_effects <- apply_latent_effects(
            n_individuals = n_individuals,
            base_capture_probabilities = base_capture_probabilities,
            latent_sd = lsd
        )
        capture_probabilities <- l_effects$capture_probabilities
        latent_effects <- l_effects$latent_effects
    } else if (het_details$latent_sd > 0) {
        l_effects <- apply_latent_effects(
            n_individuals = n_individuals,
            base_capture_probabilities = capture_probabilities,
            latent_sd = het_details$latent_sd
        )
        capture_probabilities <- l_effects$capture_probabilities
        latent_effects <- l_effects$latent_effects
    }
    return(list(
        class_assignments = class_assignments,
        latent_effects = latent_effects,
        capture_probabilities = capture_probabilities
    ))
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
