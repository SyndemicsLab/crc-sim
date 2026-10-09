################################################################################
# File: simulation_specs.R                                                     #
# Project: crcsim                                                              #
# Created Date: 2026-10-09                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-09                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Create a homogeneous capture specification
#'
#' @param p_captures Optional vector of list-specific capture probabilities.
#' @returns An S3 object with class `homogeneous_capture_spec`.
#' @export
homogeneous_capture_spec <- function(p_captures = NULL) {
    if (
        !is.null(p_captures) &&
            (!is.numeric(p_captures) ||
                length(p_captures) == 0 ||
                anyNA(p_captures) ||
                any(p_captures < 0 | p_captures > 1))
    ) {
        stop("p_captures must be NULL or probabilities between 0 and 1")
    }
    spec <- list(p_captures = p_captures)
    class(spec) <- c("homogeneous_capture_spec", "capture_spec")
    return(spec)
}

#' Create a covariate-dependent capture specification
#'
#' @param covariate_ranges Data frame describing generated covariates.
#' @param alpha List-specific logistic intercepts.
#' @param beta List-by-covariate logistic coefficient matrix.
#' @returns An S3 object with class `covariate_capture_spec`.
#' @export
covariate_capture_spec <- function(covariate_ranges, alpha, beta) {
    if (!is.data.frame(covariate_ranges) || nrow(covariate_ranges) == 0) {
        stop("covariate_ranges must be a non-empty data frame")
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
    spec <- list(
        covariate_ranges = covariate_ranges,
        alpha = alpha,
        beta = beta
    )
    class(spec) <- c("covariate_capture_spec", "capture_spec")
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
