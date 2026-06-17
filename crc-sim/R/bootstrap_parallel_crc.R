################################################################################
# File: bootstrap_parallel_crc.R                                               #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-06-17                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Return a Default List of Stepwise Selection Configurations
#' @keywords internal
#' @noRd
default_stepwise_config <- function() {
    return(list(
        f0_05 = list(direction = "forward", threshold = 0.05),
        f0_1 = list(direction = "forward", threshold = 0.1),
        b0_05 = list(direction = "backward", threshold = 0.05),
        b0_1 = list(direction = "backward", threshold = 0.1),
        fb0_05 = list(direction = "both", threshold = 0.05),
        fb0_1 = list(direction = "both", threshold = 0.1)
    ))
}

resolve_method_label <- function(method) {
    if (method == "negbin") {
        return("NB")
    }

    return("Poisson")
}

#' Run One CRC Analysis Scenario
#'
#' @keywords internal
#' @noRd
run_crc_analysis_scenario <- function(
    data_table,
    method,
    suppress,
    formula_selection,
    opts_stepwise,
    opts_tmle
) {
    captures <- names(data_table)[
        !names(data_table) %in% c("N_ID", "strata", "tmp")
    ]

    ground_truth <- extract_ground_truth(
        data_table,
        "strata",
        capture = captures
    )

    model <- simulate(
        data_table,
        "strata",
        suppress,
        method,
        formula_selection,
        opts_stepwise,
        opts_tmle,
        capture = captures
    )

    return(compare(model, ground_truth))
}

#' Run a Single CRC Method for One Config
#'
#' @keywords internal
#' @noRd
run_single_method_config <- function(
    config_value,
    data_table,
    method,
    suppress,
    formula_selection,
    opts_tmle
) {
    start_time <- proc.time()[["elapsed"]]

    output <- run_crc_analysis_scenario(
        data_table = data_table,
        method = method,
        suppress = suppress,
        formula_selection = formula_selection,
        opts_stepwise = c(config_value, verbose = FALSE),
        opts_tmle = opts_tmle
    )

    elapsed_seconds <- proc.time()[["elapsed"]] - start_time

    return(dplyr::mutate(
        output,
        time_to_convergence_sec = round(elapsed_seconds, 4)
    ))
}

#' Run Configurations for One Method
#'
#' @keywords internal
#' @noRd
run_method_across_configs <- function(
    data_table,
    method,
    config,
    suppress,
    formula_selection,
    opts_tmle
) {
    if (formula_selection == "tmle") {
        analyzed <- run_single_method_config(
            config_value = list(direction = "both", threshold = 0.05),
            data_table = data_table,
            method = method,
            suppress = suppress,
            formula_selection = formula_selection,
            opts_tmle = opts_tmle
        ) |>
            dplyr::mutate(Model = "TMLE") |>
            dplyr::mutate(Method = resolve_method_label(method))

        return(analyzed)
    }

    analyzed <- purrr::map(
        config,
        run_single_method_config,
        data_table = data_table,
        method = method,
        suppress = suppress,
        formula_selection = formula_selection,
        opts_tmle = opts_tmle
    ) |>
        format_bootstrap_analysis_output() |>
        dplyr::mutate(Method = resolve_method_label(method))

    return(analyzed)
}

#' Run Methods for One Simulated Dataset
#'
#' @keywords internal
#' @noRd
run_methods_for_one_dataset <- function(
    data_table,
    methods,
    config,
    suppress,
    formula_selection,
    opts_tmle
) {
    method_outputs <- purrr::map(
        methods,
        run_method_across_configs,
        data_table = data_table,
        config = config,
        suppress = suppress,
        formula_selection = formula_selection,
        opts_tmle = opts_tmle
    )

    return(dplyr::bind_rows(method_outputs))
}

#' Run One Single CRC Simulation
#'
#' @keywords internal
#' @noRd
run_crc_single_simulation <- function(
    p_captures,
    p_strata,
    n_individuals,
    methods,
    config,
    suppress,
    formula_selection,
    opts_tmle,
    seed
) {
    set.seed(seed)
    data_table <- create_data(n_individuals, p_captures, p_strata)

    return(run_methods_for_one_dataset(
        data_table = data_table,
        methods = methods,
        config = config,
        suppress = suppress,
        formula_selection = formula_selection,
        opts_tmle = opts_tmle
    ))
}

#' Run One Bootstrap Iteration
#'
#' @keywords internal
#' @noRd

run_crc_bootstrap_iteration <- function(
    bootstrap_index,
    n_individuals,
    p_captures,
    p_strata,
    methods,
    config,
    suppress,
    formula_selection,
    opts_tmle
) {
    invisible(bootstrap_index)
    gc() # garbage collection to free up memory between bootstraps
    data_table <- create_data(n_individuals, p_captures, p_strata)

    return(run_methods_for_one_dataset(
        data_table = data_table,
        methods = methods,
        config = config,
        suppress = suppress,
        formula_selection = formula_selection,
        opts_tmle = opts_tmle
    ))
}

#' Bind Bootstrap Output into a Single Table
#'
#' @keywords internal
#' @noRd

bind_crc_bootstrap_output <- function(output) {
    out <- dplyr::bind_rows(output, .id = "Run") |>
        dplyr::mutate(Run = as.integer(.data$Run))

    return(out)
}

#' Run CRC Simulations with Bootstrap Resampling
#'
#' @keywords internal
#' @noRd
run_crc_bootstrap_parallel <- function(
    n_bootstraps,
    p_captures,
    p_strata,
    n_individuals,
    methods,
    config,
    suppress,
    formula_selection,
    opts_tmle,
    seed
) {
    set.seed(seed)

    output <- future.apply::future_lapply(
        seq_len(n_bootstraps),
        run_crc_bootstrap_iteration,
        n_individuals = n_individuals,
        p_captures = p_captures,
        p_strata = p_strata,
        methods = methods,
        config = config,
        suppress = suppress,
        formula_selection = formula_selection,
        opts_tmle = opts_tmle,
        future.seed = TRUE
    )

    return(bind_crc_bootstrap_output(output))
}

#' Run CRC Models
#'
#' Wrapper entrypoint for package consumers. Supports either one simulation run
#' or bootstrap resampling across many runs.
#'
#' @param mode string: either "single" or "bootstrap"
#' @param n_bootstraps int: number of bootstrap runs when mode is "bootstrap"
#' @param p_captures list: list of capture probabilities
#' @param p_strata list: list of strata probabilities - for nonstratified use 1
#' @param n_individuals int: number of individuals sampled per run
#' @param methods character vector: one or both of "poisson" and "negbin"
#' @param config list: stepwise-selection configurations
#' @param suppress int: suppression cap used during simulation
#' @param formula_selection string: formula selection approach
#' @param opts_tmle list: TMLE options where \code{list_pair} is required when
#' formula_selection is \code{"tmle"}
#' @param seed int: starting seed
#'
#' @importFrom future.apply future_lapply
#' @export
run_crc <- function(
    n_individuals = 3e5,
    p_captures = c(0.45, 0.3, 0.2),
    covariate_df = 1,
    suppress = 10,
    iterations = 100,
    estimation_option = StepwiseOptions(
        direction = "both",
        threshold = 0.05,
        verbose = FALSE
    ),
    methods = c("poisson", "negbin"),
    formula_selection = "stepwise",
    seed = 2024
) {
    mode <- match.arg(mode)
    methods <- match.arg(methods, several.ok = TRUE)

    if (mode == "single") {
        return(run_crc_single_simulation(
            p_captures = p_captures,
            p_strata = p_strata,
            n_individuals = n_individuals,
            methods = methods,
            config = config,
            suppress = suppress,
            formula_selection = formula_selection,
            opts_tmle = opts_tmle,
            seed = seed
        ))
    }

    return(run_crc_bootstrap_parallel(
        n_bootstraps = n_bootstraps,
        p_captures = p_captures,
        p_strata = p_strata,
        n_individuals = n_individuals,
        methods = methods,
        config = config,
        suppress = suppress,
        formula_selection = formula_selection,
        opts_tmle = opts_tmle,
        seed = seed
    ))
}

#' Function to set up global parallel environment for future.apply
#' functionality.
#'
#' @param ncores int: The number of cores to use
#'
#' @importFrom doFuture registerDoFuture
#' @importFrom doRNG registerDoRNG
#' @importFrom future plan
#' @importFrom future multisession
#'
#' @keywords internal
#' @noRd
setup_global_parallel <- function(ncores) {
    # nolint start
    # turn off linting because we need to assign the ncores variable to the global environment for the future plan to work
    ncores <<- ncores
    # nolint end
    evalq(
        {
            doFuture::registerDoFuture()
            doRNG::registerDoRNG()
            future::plan(future::multisession(workers = ncores), gc = TRUE)
        },
        envir = .GlobalEnv
    )
    rm(ncores, envir = .GlobalEnv)
    return(NULL)
}

#' Function to format the output of the \code{analyze_wrapper} function for
#' easier comparison
#'
#' @param analyzed_dt list: list of data.tables resulting from the
#' \code{analyze_wrapper}
#'
#' @return data.table: formatted data.table with model names and method
#'
#' @keywords internal
#' @noRd
normalize_bootstrap_model_label <- function(model) {
    output <- dplyr::case_when(
        startsWith(model, "fb") ~ sub("^fb", "Both-", model),
        startsWith(model, "f") ~ sub("^f", "Forward-", model),
        startsWith(model, "b") ~ sub("^b", "Backward-", model),
        TRUE ~ model
    )

    return(output)
}

#' Format and Label Analysis Output
#'
#' @keywords internal
#' @noRd

format_bootstrap_analysis_output <- function(analyzed_dt) {
    bind_data <- dplyr::bind_rows(analyzed_dt, .id = "Model") |>
        dplyr::mutate(Model = normalize_bootstrap_model_label(.data$Model))

    return(bind_data)
}
