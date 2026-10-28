################################################################################
# File: data_mutating.R                                                        #
# Project: crcsim                                                              #
# Created Date: 2026-10-21                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-10-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' A function used to convert all integer columns to numeric. This is necessary
#' as the estimation functions require numeric data columns.
#'
#' @param df a data frame
#' @returns a data frame with all integer columns converted to numeric
#'
#' @importFrom dplyr mutate across where
#' @keywords internal
#' @export
all_int_cols_to_numeric <- function(df) {
    return(mutate(df, across(where(is.integer), as.numeric)))
}


#' Sample capture histories using sequential dependence edges.
#'
#' @param full_data Data frame containing non-capture columns and any initial
#'   capture columns.
#' @param n_pop Number of individuals to simulate.
#' @param n_capture_cols Number of capture lists.
#' @param dep_edges Directed dependence graph edges.
#' @returns A list containing the rebuilt data and final capture
#'   probabilities.
#'
#' @keywords internal
#' @noRd
update_captures_with_dependence <- function(
    full_data,
    n_pop,
    n_capture_cols,
    probs,
    dep_edges
) {
    capture_cols <- paste0("capture_", seq_len(n_capture_cols))
    res_data <- full_data[setdiff(names(full_data), capture_cols)]
    generated_captures <- matrix(0L, nrow = n_pop, ncol = n_capture_cols)

    # iterate through each capture column
    for (cap_idx in seq_len(n_capture_cols)) {
        probabilities <- probs[, cap_idx]
        in_edges <- dep_edges[dep_edges$to == cap_idx, , drop = FALSE]
        # in_edges is a dataframe with the columns `from`, `to`, and `strength`.
        # Each row represents a directed edge from `from` to `to` with a given
        # `strength`. A small example:
        #  from  to  strength
        #    1    2    0.5
        #    3    2    0.2
        if (nrow(in_edges) > 0) {
            # compute associated effect from all incoming edges. It is the sum of the effects from each incoming edge, where the effect is calculated as (capture value - 0.5) * strength and of length n_pop.
            association <- rowSums(vapply(
                seq_len(nrow(in_edges)),
                function(edge_idx) {
                    return(
                        (generated_captures[, in_edges$from[edge_idx]] - 0.5) *
                            in_edges$strength[edge_idx]
                    )
                },
                numeric(n_pop)
            ))
            probabilities <- stats::plogis(
                stats::qlogis(probabilities) + association
            )
            probs[, cap_idx] <- probabilities
        }
        generated_captures[, cap_idx] <- stats::rbinom(
            n = n_pop,
            size = 1,
            prob = probabilities
        )
    }
    for (cap_idx in seq_len(n_capture_cols)) {
        full_data[[capture_cols[[cap_idx]]]] <-
            generated_captures[, cap_idx]
    }
    full_data <- full_data[
        c(capture_cols, setdiff(names(full_data), capture_cols))
    ]
    return(full_data)
}
