################################################################################
# File: simulation_scenarios.R                                                 #
# Project: crcsim                                                              #
# Created Date: 2026-10-09                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' @title Create a Simulation Scenario
#'
#' @description This function creates a simulation scenario object that
#' encapsulates the base capture scenario, optional heterogeneity scenario,
#' and optional dependence scenario.
#'
#' @details The simulation scenario combines a base capture scenario, which can
#' be homogeneous or covariate-dependent, with optional heterogeneity and
#' dependence scenarios. The resulting object can be used to generate capture
#' histories that reflect the specified structure of the population and capture
#' process.
#'
#' @param base_scenario The base capture scenario, which can be a homogeneous
#' or covariate-dependent scenario.
#' @param heterogeneity_scenario Optional heterogeneity scenario object.
#' @param dependence_scenario Optional dependence scenario object.
#' @returns An S3 object with class `simulation_scenario`.
#'
#' @examples
#' sim_scenario <- simulation_scenario(
#'     base_scenario = homogeneous_scenario(capture_probs = rep(0.5, 6)),
#'     heterogeneity_scenario = NULL,
#'     dependence_scenario = NULL
#' )
#' This example demonstrates creating a simple simulation scenario with a
#' homogeneous base capture scenario and no additional heterogeneity or
#' dependence.
#'
#' @export
simulation_scenario <- function(
    base_scenario = homogeneous_scenario(),
    heterogeneity_scenario = NULL,
    dependence_scenario = NULL
) {
    scenario <- list(
        base_scenario = base_scenario,
        heterogeneity_scenario = heterogeneity_scenario,
        dependence_scenario = dependence_scenario
    )
    class(scenario) <- "simulation_scenario"
    return(scenario)
}

#' @title Create a homogeneous capture specification.
#'
#' @description This specifies that each capture list has a fixed probability
#' of capturing an individual, independent of any covariates or latent factors.
#'
#' @details This specification assumes each capture list is independent of the
#' others. Thus, the joint capture probability is the product of the individual
#' list probabilities. As such, the function `homogeneous_capture_spec` only
#' applies a check to the `p_captures` parameter verifying it is a numeric
#' vector of probabilities between 0 and 1 and returns it as an S3 class with
#' the `homogeneous_capture_spec` and `capture_spec` class types.
#'
#' @examples
#' hcs_one <- homogeneous_capture_spec(p_captures = rep(0.5, 6))
#' hcs_two <- homogeneous_capture_spec(p_captures = rep(0.2, 6))
#'
#' @param p_captures Optional vector of list-specific capture probabilities.
#' @returns An S3 object with class `homogeneous_capture_spec`.
#' @export
homogeneous_scenario <- function(capture_probs = NULL, n_covariates = 0) {
    if (is.null(capture_probs)) {
        capture_probs <- 0.5
    }
    validate_prob_vector(capture_probs)
    spec <- list(capture_probs = capture_probs, n_covariates = n_covariates)
    class(spec) <- c("homogeneous_scenario", "capture_scenario")
    return(spec)
}

#' @title Create a covariate-dependent capture specification.
#'
#' @description This specifies that each capture list has a probability of
#' capturing an individual depending on the covariates.
#'
#' @details This specification assumes that each capture list's probability of
#' capturing an individual is determined by a logistic model based on the
#' covariates. The lists remain conditionally independent given the covariates.
#' The function `covariate_capture_spec` checks that the inputs are valid and
#' returns them as an S3 object with the `covariate_capture_spec` and
#' `capture_spec` class types.
#'
#' \deqn{
#'  P(\text{capture}_j = 1 \mid x) =
#'      \operatorname{logit}^{-1}(\alpha_j + x^T \beta_j)
#' }
#'
#' @examples
#' ccs_one <- covariate_capture_spec(
#'     covariate_ranges = data.frame(
#'         distribution = c("uniform", "uniform"),
#'         p1 = c(0, 0),
#'         p2 = c(1, 1),
#'         dtype = c("integer", "integer")
#'     ),
#'     alpha = rep(-1, 6),
#'     beta = matrix(0.5, nrow = 6, ncol = 2)
#' )
#'
#' @param cov_dist Data frame describing generated covariates.
#' @param alpha List-specific logistic intercepts.
#' @param beta List-by-covariate logistic coefficient matrix.
#' @returns An S3 object with class `covariate_capture_spec`.
#' @export
covariate_scenario <- function(cov_dist, alpha, beta) {
    if (!is.data.frame(cov_dist) || nrow(cov_dist) == 0) {
        stop("cov_dist must be a non-empty data frame")
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
    if (ncol(beta) != nrow(cov_dist)) {
        stop("beta must have one column for each covariate")
    }
    spec <- list(
        cov_dists = cov_dists,
        alpha = alpha,
        beta = beta
    )
    class(spec) <- c("covariate_scenario", "capture_scenario")
    return(spec)
}

#' Create a continuous latent heterogeneity specification
#'
#' @param sd Standard deviation of the shared latent propensity.
#' @returns An S3 object with class `continuous_heterogeneity_spec`.
#' @export
continuous_heterogeneity_spec <- function(sd) {
    if (!is.numeric(sd) || length(sd) != 1 || is.na(sd) || sd < 0) {
        stop("sd must be a non-negative numeric scalar")
    }
    spec <- list(sd = sd)
    class(spec) <- "continuous_heterogeneity_spec"
    return(spec)
}

#' Create a latent-class heterogeneity specification
#'
#' @param probabilities Latent-class assignment probabilities.
#' @param capture_probabilities Class-by-capture probability matrix.
#' @returns An S3 object with class `latent_class_spec`.
#' @export
latent_class_spec <- function(probabilities, capture_probabilities) {
    if (
        !is.numeric(probabilities) ||
            length(probabilities) < 2 ||
            anyNA(probabilities) ||
            any(probabilities < 0) ||
            sum(probabilities) <= 0
    ) {
        stop("probabilities must be non-negative and sum to a positive value")
    }
    if (
        !is.numeric(capture_probabilities) ||
            !is.matrix(capture_probabilities) ||
            anyNA(capture_probabilities) ||
            any(capture_probabilities < 0 | capture_probabilities > 1)
    ) {
        stop("capture_probabilities must be a probability matrix")
    }
    if (nrow(capture_probabilities) != length(probabilities)) {
        stop("capture_probabilities needs one row per latent class")
    }
    spec <- list(
        probabilities = probabilities / sum(probabilities),
        capture_probabilities = capture_probabilities
    )
    class(spec) <- "latent_class_spec"
    return(spec)
}

#' Create a pairwise dependence specification
#'
#' @param pair Ordered pair of dependent capture-list indices.
#' @param strength Conditional-logit dependence strength.
#' @returns An S3 object with class `pairwise_dependence_spec`.
#' @export
pairwise_dependence_spec <- function(pair, strength = 0) {
    validate_spec_dependence_pair(pair)
    if (!is.numeric(strength) || length(strength) != 1 || is.na(strength)) {
        stop("strength must be a numeric scalar")
    }
    spec <- list(pair = pair, strength = strength)
    class(spec) <- "pairwise_dependence_spec"
    return(spec)
}

#' Create a graph dependence specification
#'
#' @param edges Data frame with `from`, `to`, and `strength` columns.
#' @returns An S3 object with class `graph_dependence_spec`.
#' @export
graph_dependence_spec <- function(edges) {
    validate_spec_dependence_edges(edges)
    spec <- list(edges = edges)
    class(spec) <- "graph_dependence_spec"
    return(spec)
}

validate_spec_dependence_pair <- function(pair) {
    if (
        !is.numeric(pair) ||
            length(pair) != 2 ||
            anyNA(pair) ||
            any(pair %% 1 != 0) ||
            any(pair < 1) ||
            pair[1] >= pair[2]
    ) {
        stop("pair must be an ordered pair of distinct list indices")
    }
    return(invisible(NULL))
}

validate_spec_dependence_edges <- function(edges) {
    if (
        !is.data.frame(edges) ||
            !all(c("from", "to", "strength") %in% names(edges)) ||
            nrow(edges) < 1
    ) {
        stop("edges must contain at least one from, to, and strength row")
    }
    if (
        !is.numeric(edges$from) ||
            !is.numeric(edges$to) ||
            !is.numeric(edges$strength) ||
            anyNA(edges) ||
            any(edges$from %% 1 != 0) ||
            any(edges$to %% 1 != 0) ||
            any(edges$from < 1) ||
            any(edges$to < 1) ||
            any(edges$from >= edges$to)
    ) {
        stop("edges must contain ordered numeric graph edges")
    }
    return(invisible(NULL))
}
