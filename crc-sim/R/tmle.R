################################################################################
# File: tmle.R                                                                 #
# Project: crc-sim                                                             #
# Created Date: 2026-05-13                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-13                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Validate TMLE Options
#'
#' @keywords internal
#' @noRd

validate_tmle_opts <- function(opts_tmle, binary_variables) {
    if (is.null(opts_tmle$list_pair)) {
        stop("For formula_selection = 'tmle', opts_tmle$list_pair is required")
    }

    if (
        !is.character(opts_tmle$list_pair) || length(opts_tmle$list_pair) != 2
    ) {
        stop("opts_tmle$list_pair must be a character vector of length 2")
    }

    if (opts_tmle$list_pair[1] == opts_tmle$list_pair[2]) {
        stop("opts_tmle$list_pair must contain two distinct capture columns")
    }

    if (!all(opts_tmle$list_pair %in% binary_variables)) {
        stop("opts_tmle$list_pair must be columns in binary_variables")
    }

    expansion_mode <- opts_tmle$expansion_mode
    if (!(expansion_mode %in% c("exact", "sample"))) {
        stop("opts_tmle$expansion_mode must be 'exact' or 'sample'")
    }

    if (
        expansion_mode == "sample" &&
            (is.null(opts_tmle$sample_size) || opts_tmle$sample_size < 2)
    ) {
        stop(
            "opts_tmle$sample_size must be >= 2 when expansion_mode = 'sample'"
        )
    }

    return(NULL)
}

#' Expand Contingency Table for TMLE Input
#'
#' @keywords internal
#' @noRd

expand_tmle_input <- function(data, freq_column, opts_tmle) {
    freq <- data[[freq_column]]

    if (any(freq < 0) || anyNA(freq)) {
        stop("Frequency column must contain non-negative, non-missing values")
    }

    if (!all(freq %% 1 == 0)) {
        stop("Frequency column must contain integer counts")
    }

    observed_n <- sum(freq)
    expansion_mode <- opts_tmle$expansion_mode

    if (expansion_mode == "exact") {
        if (observed_n > opts_tmle$warn_expanded_rows) {
            warning(paste(
                "TMLE exact expansion generated",
                observed_n,
                "rows; this may be memory intensive"
            ))
        }

        expanded <- data |>
            dplyr::slice(
                rep(dplyr::row_number(), times = freq)
            ) |>
            tibble::as_tibble()

        return(list(expanded = expanded, observed_n = observed_n))
    }

    row_count <- nrow(data)
    sampled_index <- sample(
        seq_len(row_count),
        size = opts_tmle$sample_size,
        replace = TRUE,
        prob = freq / observed_n
    )

    expanded <- data |>
        dplyr::slice(sampled_index) |>
        tibble::as_tibble()

    return(list(expanded = expanded, observed_n = observed_n))
}

#' Resolve TMLE Result Row
#'
#' @keywords internal
#' @noRd

resolve_tmle_result_row <- function(psin_estimate, opts_tmle, j, k) {
    pair_label <- paste(min(j, k), max(j, k), sep = ",")
    result <- psin_estimate$result

    preferred_model <- opts_tmle$funcname[1]
    preferred_method <- opts_tmle$estimator

    pair_subset <- dplyr::filter(result, listpair == pair_label)
    if (nrow(pair_subset) == 0) {
        stop("TMLE returned no result for selected list_pair")
    }

    model_subset <- dplyr::filter(pair_subset, model == preferred_model)
    if (nrow(model_subset) == 0) {
        model_subset <- pair_subset
    }

    method_subset <- dplyr::filter(model_subset, method == preferred_method)
    if (nrow(method_subset) == 0) {
        stop("TMLE returned no result for selected estimator method")
    }

    return(dplyr::slice(method_subset, 1))
}

#' Run TMLE for CRC
#'
#' @keywords internal
#' @noRd

run_crc_tmle <- function(data, freq_column, binary_variables, opts_tmle) {
    validate_tmle_opts(opts_tmle, binary_variables)

    expanded_obj <- expand_tmle_input(data, freq_column, opts_tmle)
    expanded_data <- expanded_obj$expanded
    observed_n <- expanded_obj$observed_n

    covariate_columns <- names(expanded_data) |>
        setdiff(c(freq_column, binary_variables))

    tmle_data <- expanded_data |>
        dplyr::mutate(
            dplyr::across(
                dplyr::all_of(binary_variables),
                as.numeric
            )
        ) |>
        dplyr::select(
            dplyr::all_of(c(binary_variables, covariate_columns))
        )

    j <- match(opts_tmle$list_pair[1], binary_variables)
    k <- match(opts_tmle$list_pair[2], binary_variables)

    psin_estimate <- drpop::popsize(
        data = tmle_data,
        K = length(binary_variables),
        j = j,
        k = k,
        margin = opts_tmle$margin,
        nfolds = opts_tmle$nfolds,
        funcname = opts_tmle$funcname,
        TMLE = TRUE,
        PLUGIN = TRUE
    )

    result_row <- resolve_tmle_result_row(psin_estimate, opts_tmle, j, k)

    estimate_hidden <- as.numeric(result_row$n) - observed_n
    lower_hidden <- as.numeric(result_row$cin.l) - observed_n
    upper_hidden <- as.numeric(result_row$cin.u) - observed_n

    output <- list(
        model = "tmle",
        formula = paste("pair:", paste(opts_tmle$list_pair, collapse = ",")),
        summary = psin_estimate,
        estimate = round(max(estimate_hidden, 0), 2),
        lower_ci = round(max(lower_hidden, 0), 2),
        upper_ci = round(max(upper_hidden, 0), 2),
        AIC = NA_real_,
        observed_n = observed_n,
        tmle_method = as.character(result_row$method),
        tmle_model = as.character(result_row$model),
        list_pair = as.character(result_row$listpair)
    )

    return(output)
}
