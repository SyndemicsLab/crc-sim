################################################################################
# File: petersen_estimate.R                                                    #
# Project: crc-sim                                                             #
# Created Date: 2026-05-18                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-20                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Compute Petersen estimates for all pairs of capture indicators in the data.
#'
#' This function takes an estimated q1, q2, q12 for each pair of capture
#' indicators and computes Petersen estimates and standard errors for each
#' pair. It returns a list of results for each pair.
#' @param q1 Vector of predicted capture probabilities for list j.
#' @param q2 Vector of predicted capture probabilities for list k.
#' @param q12 Vector of predicted joint capture probabilities for lists j and k.
#' @param n_obs Number of observations in the test set.
#' @return A list of Petersen estimates and standard errors for each pair of
#' capture indicators.
#' @keywords internal
#' @export
petersen_pair_estimate <- function(
    q1,
    q2,
    q12,
    n_obs
) {
    if (q12 == 0) {
        warning(
            paste(
                "No overlap between capture indicators",
                "- cannot compute Petersen estimate."
            )
        )
        return(list(
            n_hat = NA,
            sd_n_hat = NA,
            psi_inverse_hat = NA,
            se_psi_inverse_hat = NA
        ))
    }

    psi_inverse_hat <- max(q1 * q2 / q12, 1)
    se_psi_inverse_hat <- sqrt(
        q1 *
            q2 *
            pmax(q1 * q2 - q12, 0) *
            (1 - q12) /
            q12^3 /
            n_obs
    )

    sd_n_hat <- sqrt(
        n_obs^2 *
            se_psi_inverse_hat^2 +
            n_obs * psi_inverse_hat * (psi_inverse_hat - 1)
    )

    return(list(
        n_hat = round(n_obs * psi_inverse_hat),
        sd_n_hat = sd_n_hat,
        psi_inverse_hat = psi_inverse_hat,
        se_psi_inverse_hat = se_psi_inverse_hat
    ))
}

#' Compute Petersen estimates for all pairs of capture indicators in the data.
#'
#' This function takes a data frame with binary capture indicators as columns
#' and computes Petersen estimates and standard errors for each pair of capture
#' indicators. It returns a list of results for each pair.
#' @param data A data frame containing binary capture indicators as columns.
#' @return A list of Petersen estimates and standard errors for each pair of
#' capture indicators.
#' @importFrom utils combn
#' @export
petersen_estimate <- function(data) {
    if (!is.data.frame(data)) {
        stop("`data` must be a data frame.")
    }

    if (ncol(data) < 2) {
        stop("`data` must contain at least two capture columns.")
    }
    if (
        lapply(data, function(col) !all(col %in% c(0, 1))) |>
            unlist() |>
            any()
    ) {
        stop(paste0(
            "All columns in `data` must contain only 0 or 1",
            "capture indicators."
        ))
    }
    column_names <- names(data)
    pair_indices <- combn(column_names, 2, simplify = FALSE)

    pair_results <- lapply(pair_indices, function(pair) {
        left <- pair[[1]]
        right <- pair[[2]]

        x_left <- data[[left]]
        x_right <- data[[right]]
        n_obs <- nrow(data)

        q1 <- mean(x_left)
        q2 <- mean(x_right)
        q12 <- mean(x_left * x_right)
        res <- petersen_pair_estimate(
            q1 = q1,
            q2 = q2,
            q12 = q12,
            n_obs = n_obs
        )
        return(res)
    })

    names(pair_results) <- vapply(
        pair_indices,
        paste,
        collapse = ",",
        character(1)
    )
    return(pair_results)
}
