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
#' @param n_captures int: optional number of binary capture columns to simulate.
#'   In covariate mode, this is inferred from \\code{alpha} when omitted.
#' @param p_captures numeric vector: optional list-specific capture
#'   probabilities. Defaults to 0.5 for each list.
#' @param latent_class_probabilities numeric vector: optional probabilities for
#'   unobserved latent capture classes.
#' @param latent_class_p_captures numeric matrix: optional class-by-list matrix
#'   of capture probabilities. Rows correspond to latent classes.
#' @param covariate_ranges data.frame: optional covariate specification passed
#'   to \code{create_data}.
#' @param alpha numeric vector: optional list-specific logistic intercepts.
#' @param beta numeric matrix: optional list-by-covariate logistic coefficient
#'   matrix.
#' @param latent_sd numeric scalar: standard deviation of a shared latent
#'   capture propensity on the logit scale. Defaults to zero.
#' @param dependence_pair integer vector of length two: optional ordered pair
#'   of capture lists with direct conditional dependence.
#' @param dependence_strength numeric scalar: conditional-logit association
#'   strength for \\code{dependence_pair}. Defaults to zero.
#' @param dependence_edges data.frame: optional directed graph edges with
#'   integer columns \\code{from}, \\code{to}, and numeric \\code{strength}.
#'   Edges must point from a lower-numbered list to a higher-numbered list.
#' @param seed numeric scalar: optional random seed for reproducible output.
#' @details
#' The returned object separates the complete generated population from the
#' observed union and the all-zero population component. The all-zero component
#' is retained as simulation truth and is not available from an observed-only
#' data frame. Fixed-probability simulations record their normalized list
#' probabilities in \\code{truth$p_captures}; covariate-dependent simulations
#' record one true probability per individual and capture list in
#' \\code{truth$capture_probabilities}. A positive \\code{latent_sd} adds a
#' shared unobserved propensity, while \\code{dependence_pair} changes the
#' second list's conditional log odds based on the first list's capture.
#' Latent-class mode assigns each individual to an unobserved class and samples
#' lists independently conditional on that class.
#' @returns A list with `full_data`, `observed_data`, `unobserved_data`,
#'   `observed_frequency_table`, `truth`, and `metadata` components.
#'
#' @export
simulate_data <- function(
    n_individuals,
    n_captures = NULL,
    p_captures = NULL,
    latent_class_probabilities = NULL,
    latent_class_p_captures = NULL,
    covariate_ranges = NULL,
    alpha = NULL,
    beta = NULL,
    latent_sd = 0,
    dependence_pair = NULL,
    dependence_strength = 0,
    dependence_edges = NULL,
    seed = NULL
) {
    validate_simulation_seed(seed)
    validate_dependence_args(
        n_captures = n_captures,
        latent_sd = latent_sd,
        dependence_pair = dependence_pair,
        dependence_strength = dependence_strength,
        dependence_edges = dependence_edges
    )
    validate_latent_class_args(
        n_captures = n_captures,
        p_captures = p_captures,
        alpha = alpha,
        beta = beta,
        latent_class_probabilities = latent_class_probabilities,
        latent_class_p_captures = latent_class_p_captures
    )

    if (!is.null(seed)) {
        set.seed(seed)
    }

    covariate_mode <- !is.null(alpha) || !is.null(beta)
    latent_class_mode <- !is.null(latent_class_probabilities)
    if (latent_class_mode && is.null(n_captures)) {
        n_captures <- ncol(latent_class_p_captures)
    }
    if (latent_class_mode) {
        latent_class_probabilities <- latent_class_probabilities /
            sum(latent_class_probabilities)
    }
    class_assignments <- NULL
    latent_effects <- rep(0, n_individuals)
    mechanism_override <- latent_sd > 0 ||
        !is.null(dependence_pair) ||
        latent_class_mode ||
        !is.null(dependence_edges)
    graph_mode <- !is.null(dependence_edges)
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
    if (!covariate_mode) {
        if (is.null(n_captures)) {
            stop("n_captures is required without alpha and beta")
        }
        full_data <- create_data(
            n_individuals = n_individuals,
            n_captures = n_captures,
            p_captures = p_captures,
            covariate_ranges = covariate_ranges
        )
        capture_probabilities <- NULL
        p_captures <- if (is.null(p_captures)) {
            rep(0.5, n_captures)
        } else {
            as.numeric(p_captures)
        }
        scenario <- "homogeneous conditional independence"
        base_capture_probabilities <- matrix(
            rep(p_captures, each = n_individuals),
            nrow = n_individuals,
            ncol = n_captures
        )
    } else {
        validate_covariate_capture_model(
            n_captures = n_captures,
            p_captures = p_captures,
            covariate_ranges = covariate_ranges,
            alpha = alpha,
            beta = beta
        )
        n_captures <- length(alpha)
        covariate_data <- create_data(
            n_individuals = n_individuals,
            n_captures = 1,
            p_captures = 0,
            covariate_ranges = covariate_ranges
        )
        full_data <- covariate_data[,
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
        full_data <- full_data[,
            c(capture_columns, covariate_columns),
            drop = FALSE
        ]
        scenario <- "covariate-dependent conditional independence"
    }

    if (latent_class_mode) {
        n_captures <- ncol(latent_class_p_captures)
        class_assignments <- sample(
            seq_along(latent_class_probabilities),
            size = n_individuals,
            replace = TRUE,
            prob = latent_class_probabilities
        )
        base_capture_probabilities <- latent_class_p_captures[
            class_assignments,
            ,
            drop = FALSE
        ]
    }

    capture_columns <- paste0("capture_", seq_len(n_captures))
    if (mechanism_override) {
        if (latent_sd > 0) {
            latent_effects <- stats::rnorm(n_individuals, sd = latent_sd)
            capture_probabilities <- stats::plogis(
                stats::qlogis(base_capture_probabilities) + latent_effects
            )
        } else {
            capture_probabilities <- base_capture_probabilities
        }
        full_data <- full_data[,
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
                        (generated_captures[,
                            incoming_edges$from[edge_index]
                        ] -
                            0.5) *
                            incoming_edges$strength[edge_index]
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
        full_data <- full_data[,
            c(capture_columns, setdiff(names(full_data), capture_columns)),
            drop = FALSE
        ]
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
    }
    observed_data <- extract_captured_data(full_data, capture_columns)
    unobserved_data <- extract_uncaptured_data(full_data, capture_columns)

    result <- list(
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
    )
    class(result) <- c("crcsim_simulation", "list")
    return(result)
}

validate_dependence_args <- function(
    n_captures,
    latent_sd,
    dependence_pair,
    dependence_strength,
    dependence_edges
) {
    if (
        !is.numeric(latent_sd) ||
            length(latent_sd) != 1 ||
            is.na(latent_sd) ||
            latent_sd < 0
    ) {
        stop("latent_sd must be a non-negative numeric scalar")
    }
    if (
        !is.numeric(dependence_strength) ||
            length(dependence_strength) != 1 ||
            is.na(dependence_strength)
    ) {
        stop("dependence_strength must be a numeric scalar")
    }
    if (!is.null(dependence_pair) && !is.null(dependence_edges)) {
        stop("supply either dependence_pair or dependence_edges, not both")
    }
    if (is.null(dependence_pair) && is.null(dependence_edges)) {
        return(invisible(NULL))
    }
    if (!is.null(dependence_edges)) {
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
        if (
            !is.null(n_captures) &&
                any(dependence_edges$to > n_captures)
        ) {
            stop("dependence_edges must refer to existing capture lists")
        }
        return(invisible(NULL))
    }
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
    if (!is.null(n_captures) && dependence_pair[2] > n_captures) {
        stop("dependence_pair must refer to existing capture lists")
    }
    return(invisible(NULL))
}

validate_latent_class_args <- function(
    n_captures,
    p_captures,
    alpha,
    beta,
    latent_class_probabilities,
    latent_class_p_captures
) {
    if (
        is.null(latent_class_probabilities) &&
            is.null(latent_class_p_captures)
    ) {
        return(invisible(NULL))
    }
    if (
        is.null(latent_class_probabilities) ||
            is.null(latent_class_p_captures)
    ) {
        stop(
            "latent_class_probabilities and latent_class_p_captures must ",
            "be supplied together"
        )
    }
    if (!is.null(alpha) || !is.null(beta)) {
        stop("latent classes cannot be combined with alpha and beta")
    }
    if (!is.null(p_captures)) {
        stop("p_captures cannot be supplied with latent classes")
    }
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
    if (nrow(latent_class_p_captures) != length(latent_class_probabilities)) {
        stop("latent_class_p_captures needs one row per latent class")
    }
    if (!is.null(n_captures) && ncol(latent_class_p_captures) != n_captures) {
        stop("latent_class_p_captures needs one column per capture list")
    }
    return(invisible(NULL))
}

validate_simulation_seed <- function(seed) {
    if (
        !is.null(seed) &&
            (!is.numeric(seed) ||
                length(seed) != 1 ||
                is.na(seed) ||
                seed %% 1 != 0)
    ) {
        stop("seed must be NULL or an integer scalar")
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
    if (!is.null(n_captures)) {
        if (
            !is.numeric(n_captures) ||
                length(n_captures) != 1 ||
                n_captures < 1 ||
                n_captures %% 1 != 0
        ) {
            stop("n_captures must be a positive integer")
        }
        if (n_captures != length(alpha)) {
            stop("n_captures must equal the length of alpha")
        }
    }
    return(invisible(NULL))
}
