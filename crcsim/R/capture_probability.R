################################################################################
# File: capture_probability.R                                                  #
# Project: crcsim                                                              #
# Created Date: 2026-08-26                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Estimate capture probability with or without covariates.
#'
#' @param data The dataset containing capture columns and optionally covariates.
#' @param n_lists The number of capture lists.
#' @param method The estimation method to be used (default is "doubly_robust").
#' @param func The function used for nuisance parameter estimation (default is
#' "logit").
#' @param nfolds The number of cross-validation folds (default is 5).
#' @param margin The margin parameter for estimation (default is 0.005).
#' @param seed The random seed for reproducibility (default is NULL).
#' @param ... Additional arguments passed to the nuisance estimation function.
#' @return A list containing the estimation results and the total population
#' size.
#'
#' @export
estimate_capture_prob <- function(
    data,
    n_lists,
    method = "doubly_robust",
    func = "logit",
    nfolds = 5,
    margin = 0.005,
    seed = NULL,
    ...
) {
    # validate format first: check column names/order and no NA values
    stopifnot(validate_binary_cols(data, n_lists))

    n_covariates <- ncol(data) - n_lists
    n_obs <- nrow(data)

    if (nfolds > 1 && nfolds > n_obs / 50) {
        nfolds <- pmax(floor(n_obs / 50), 1)
        warning(paste0(
            "nfolds is reduced to ",
            nfolds,
            " to have sufficient test data.\n"
        ))
    }

    if (n_obs < 1) {
        stop("Data must have at least one row.")
    }

    if (n_covariates == 0) {
        return(estimate_no_covariates(data, n_lists))
    }
    return(estimate_with_covariates(
        data,
        n_lists,
        method,
        func,
        nfolds,
        margin,
        seed
    ))
}


#' Estimate population size without covariates.
#'
#' @param data The dataset containing capture columns.
#' @param n_lists The number of capture lists.
#' @return A list containing the estimation results and the total population
#' size.
#'
#' @importFrom dplyr mutate select
#' @importFrom tidyr separate_wider_delim
#' @importFrom purrr map2_dbl
#'
#' @export
estimate_no_covariates <- function(data, n_lists) {
    n_obs <- nrow(data)
    # revalidate that all columns are binary 0/1 values
    stopifnot(validate_binary_cols(data, n_lists))

    listpairs <- build_listpairs(n_lists)

    estimates <- data.frame(listpair = listpairs) |>
        separate_wider_delim(
            listpair,
            delim = ",",
            names = c("j", "k")
        ) |>
        mutate(
            listpair = paste0(j, ",", k),
            q1 = map2_dbl(j, k, ~ mean(data[[as.integer(.x)]])),
            q2 = map2_dbl(j, k, ~ mean(data[[as.integer(.y)]])),
            q12 = map2_dbl(
                j,
                k,
                ~ mean(data[[as.integer(.x)]] * data[[as.integer(.y)]])
            ),
            psi_inv = pmax(q1 * q2 / q12, 1),
            sigma = sqrt(
                q1 * q2 * pmax(q1 * q2 - q12, 0) *
                    (1 - q12) / q12^3 / n_obs
            )
        ) |>
        select(listpair, psi_inv, sigma) |>
        mutate(
            sigma = sqrt(n_obs) * sigma,
            n = round(n_obs * psi_inv),
            sigma_n = sqrt(
                n_obs^2 * sigma^2 +
                    n_obs * psi_inv * (psi_inv - 1)
            ),
            ci_l = round(pmax(
                n_obs *
                    psi_inv -
                    1.96 * sqrt(
                        n_obs^2 * sigma^2 +
                            n_obs * psi_inv * (psi_inv - 1)
                    ),
                n_obs
            )),
            ci_u = round(
                n_obs *
                    psi_inv +
                    1.96 * sqrt(
                        n_obs^2 * sigma^2 +
                            n_obs * psi_inv * (psi_inv - 1)
                    )
            )
        ) |>
        select(listpair, n, sigma_n, ci_l, ci_u)
    return(estimates)
}

#' Estimate population size with covariates using cross-validated nuisance
#' parameter estimation.
#'
#' @param method The estimation method to be used.
#' @param data The dataset containing capture columns and covariates.
#' @param n_lists The number of capture lists.
#' @param funcname The function used for nuisance parameter estimation.
#' @param seed The random seed for reproducibility (default is NULL).
#' @return A list containing the estimation results and the total population
#' size.
#'
#' @export
estimate_with_covariates <- function(
    data,
    n_lists,
    method,
    func,
    nfolds,
    margin,
    seed = NULL
) {
    ############################################################################
    # Crossfold Setup
    ############################################################################
    if (!is.null(seed)) {
        stopifnot(length(seed) == 1, is.numeric(seed), !is.na(seed))
        set.seed(seed)
    }
    permutset <- sample(seq_len(nrow(data)), nrow(data), replace = FALSE)

    # Creates list of list(train, test, sbset)
    folds <- lapply(
        seq_len(nfolds),
        build_folds,
        data = data,
        nfolds = nfolds,
        permutation = permutset
    )

    ############################################################################
    # The Estimation Process
    ############################################################################

    listpairs <- build_listpairs(n_lists)
    list1_vec <- 1:(n_lists - 1)
    list2_vec <- 2:n_lists

    n_obs <- nrow(data)

    summaries <- purrr::map_dfr(
        listpairs,
        summarize_pair,
        folds = folds,
        method = method,
        margin = margin,
        n_lists = n_lists,
        n = n_obs,
        func = func
    )

    estimates <- summaries |>
        dplyr::mutate(
            sigma = sqrt(n_obs * var),
            n = round(n_obs * psi_inverse),
            sigma_n = sqrt(
                n_obs^2 * var +
                    n_obs * psi_inverse * (psi_inverse - 1)
            ),
            ci_l = round(pmax(
                n_obs *
                    psi_inverse -
                    1.96 *
                        sqrt(
                            n_obs^2 * var +
                                n_obs * psi_inverse * (psi_inverse - 1)
                        ),
                n_obs
            )),
            ci_u = round(
                n_obs *
                    psi_inverse +
                    1.96 *
                        sqrt(
                            n_obs^2 * var +
                                n_obs * psi_inverse * (psi_inverse - 1)
                        )
            )
        ) |>
        select(listpair, n, sigma_n, ci_l, ci_u)

    return(estimates)
}

#' Summarize the results for a specific pair of capture lists across all
#' cross-validation folds.
#'
#' @param listpair A character string representing the pair of capture lists
#' (e.g., "1,2").
#' @param folds A list of cross-validation folds, each containing training and
#' test sets.
#' @param ... Additional arguments passed to the `run_fold` function.
#' @return A tibble summarizing the capture probability and variance for the
#' specified pair of capture lists.
#'
#' @importFrom purrr map map_dbl compact
#' @importFrom tibble tibble
#' @keywords internal
summarize_pair <- function(listpair, folds, ...) {
    parts <- as.integer(strsplit(listpair, ",", fixed = TRUE)[[1]])

    fold_results <- map(
        folds,
        run_fold,
        j = parts[[1]],
        k = parts[[2]],
        ...
    ) |>
        compact()

    summary_table <- tibble(
        listpair = listpair,
        psi_inverse = mean(
            map_dbl(fold_results, "capture_probability")
        ),
        var = mean(
            map_dbl(fold_results, "variance")
        )
    )
    return(summary_table)
}

#' Estimate the capture probability and its variance using the specified method.
#'
#' @param method A character string specifying the estimation method ("tmle" or
#' "doubly_robust").
#' @param nuisance_functions A list containing the nuisance functions q_1, q_2,
#' and q_12.
#' @param y_j A numeric vector representing the observations for capture list j.
#' @param y_k A numeric vector representing the observations for capture list k.
#' @param margin A numeric value specifying the margin parameter.
#' @param n_lists An integer specifying the number of capture lists.
#' @param n An integer specifying the sample size.
#' @return A list containing the estimated capture probability and its variance.
#'
#' @keywords internal
estimate_psi <- function(
    method,
    nuisance_functions,
    y_j,
    y_k,
    margin,
    n_lists,
    n
) {
    q_1 <- nuisance_functions[["q_1"]]
    q_2 <- nuisance_functions[["q_2"]]
    q_12 <- nuisance_functions[["q_12"]]
    if (method == "tmle") {
        updated_nuisances <- tmle_nuisance(
            q_1,
            q_2,
            q_12,
            y_j,
            y_k,
            y_j * y_k,
            iterations = 100,
            margin = margin,
            n_lists = n_lists
        )
        q_1 <- updated_nuisances[["q_1"]]
        q_2 <- updated_nuisances[["q_2"]]
        q_12 <- updated_nuisances[["q_12"]]
    }

    q_1 <- pmin(pmax(q_12, q_1), 1)
    q_2 <- pmax(q_12 / q_1, pmin(q_2, 1 + q_12 - q_1, 1))

    phi_hat <- influence(y_j, y_k, q_1, q_2, q_12)
    variance <- var(phi_hat, na.rm = TRUE) / n

    if (method == "doubly_robust") {
        capture_probability <- doubly_robust_estimator(
            q_1,
            q_2,
            q_12,
            phi_hat
        )
    } else {
        capture_probability <- plugin_estimator(q_1, q_2, q_12)
    }

    return(list(
        capture_probability = capture_probability,
        variance = variance
    ))
}
