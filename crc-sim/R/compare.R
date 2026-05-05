################################################################################
# File: compare.R                                                              #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-05                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Build One Comparison Row
#'
#' @param sim_result list: one simulation result with AIC and estimate fields
#' @param ground_count integer: uncaptured ground-truth count for one group
#' @param group_name character: label for the corresponding group
#' @returns tibble with one comparison row
#' @keywords internal
#' @noRd
build_compare_row <- function(sim_result, ground_count, group_name) {
    return(tibble::tibble(
        aic = sim_result$AIC,
        estimate = sim_result$estimate,
        lower_ci = sim_result$lower_ci,
        upper_ci = sim_result$upper_ci,
        ground = ground_count,
        group = group_name
    ))
}

#' Build Grouped Comparison Output
#'
#' @param sim list: grouped simulation outputs from \code{simulate}
#' @param ground_truth list: grouped uncaptured ground-truth values
#' @returns tibble with one row per group
#' @keywords internal
#' @noRd
build_grouped_compare_output <- function(sim, ground_truth) {
    group_names <- sort(names(ground_truth))
    compare_rows <- purrr::map(
        seq_along(sim),
        function(i) {
            return(build_compare_row(
                sim_result = sim[[i]],
                ground_count = as.integer(ground_truth[[group_names[i]]]),
                group_name = group_names[i]
            ))
        }
    )

    return(dplyr::bind_rows(compare_rows))
}

#' Build Base Comparison Output
#'
#' @param sim list: ungrouped simulation output from \code{simulate}
#' @param ground_truth list: ungrouped uncaptured ground-truth value
#' @returns tibble with a single "base" comparison row
#' @keywords internal
#' @noRd
build_base_compare_output <- function(sim, ground_truth) {
    return(build_compare_row(
        sim_result = sim,
        ground_count = as.integer(ground_truth[1]),
        group_name = "base"
    ))
}

#' Compare Simulation Results Against Ground Truth
#'
#' @param sim list: simulations produced by \code{simulate}
#' @param ground_truth list: ground truth values produced by
#' \code{extract_ground_truth}
#' @returns tibble with AIC, estimate, confidence intervals, and ground truth
#' for each group
#' @export

compare <- function(sim, ground_truth) {
    if (length(ground_truth) > 1) {
        return(build_grouped_compare_output(sim, ground_truth))
    }

    return(build_base_compare_output(sim, ground_truth))
}
