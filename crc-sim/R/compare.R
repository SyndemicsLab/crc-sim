################################################################################
# File: compare.R                                                              #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-02-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Function to compare between simulation results and ground truth
#'
#' @param sim list of simulations done by \code{simulate}
#' @param ground_truth list of ground truths by \code{extract_ground_truth}
#'
#' @import data.table
#' @keywords internal

compare <- function(sim, ground_truth) {
    compare_list <- c()
    if (length(ground_truth) > 1) {
        for (i in seq_along(sim)) {
            aic <- sim[[i]]$AIC
            est <- sim[[i]]$estimate
            est_lci <- sim[[i]]$lower_ci
            est_uci <- sim[[i]]$upper_ci

            group <- sort(names(ground_truth))[i]
            gt <- as.integer(ground_truth[[group]])

            out_list <- list(
                aic = aic,
                estimate = est,
                lower_ci = est_lci,
                upper_ci = est_uci,
                ground = gt,
                group = group
            )
            compare_list[[i]] <- out_list
        }
        out <- rbindlist(compare_list)
    } else {
        aic <- sim$AIC
        est <- sim$estimate
        est_lci <- sim$lower_ci
        est_uci <- sim$upper_ci
        group <- "base"
        gt <- as.integer(ground_truth[1])

        out <- data.table(
            aic = aic,
            estimate = est,
            lower_ci = est_lci,
            upper_ci = est_uci,
            ground = gt,
            group = group
        )
    }
    return(out)
}
