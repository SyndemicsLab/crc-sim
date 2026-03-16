################################################################################
# File: run_crc.R                                                              #
# Project: crc-sim                                                             #
# Created Date: 2026-02-17                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-02-23                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Function to run CRC sim - returns a list of data.tables resulting from
#' \code{analyze}
#'
#' @param nboot int: number of bootstraps
#' @param ncores int: number of cores to use
#' @param p_captures list: list of capture probabilities
#' @param p_strata list: list of strata probabilities - for nonstratified just
#' use '1'
#' @param seed int: starting seed
#'
#' @importFrom future.apply future_lapply
#' @import data.table
#'
#' @export
run_crc <- function(
    n_bootstraps,
    p_captures,
    p_strata,
    n_individuals = 3e5,
    config = list(
        f0_05 = list(direction = "forward", threshold = 0.05),
        f0_1 = list(direction = "forward", threshold = 0.1),
        b0_05 = list(direction = "backward", threshold = 0.05),
        b0_1 = list(direction = "backward", threshold = 0.1),
        fb0_05 = list(direction = "both", threshold = 0.05),
        fb0_1 = list(direction = "both", threshold = 0.1)
    ),
    seed = 2024
) {
    set.seed(seed)
    output <- future.apply::future_lapply(
        1:n_bootstraps,
        function(x) {
            gc() # garbage collection to free up memory between bootstraps
            return(run_poisson_and_negbin(
                n_individuals,
                p_captures,
                p_strata,
                config
            ))
        },
        future.seed = TRUE
    )
    return(rbindlist(output, idcol = "Run"))
}


#' Function to set up global parallel environment for future.apply
#' functionality.
#'
#' @param ncores int: The number of cores to use
#'
#' @importFrom doFuture registerDoFuture
#' @importFrom doRNG registerDoRNG
#' @importFrom future plan_multisession
#'
#' @keywords internal
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

#' Wrapper function to run the analyze function with the same parameters.
#'
#' @param x list: list of parameters for stepwise selection - direction and
#' threshold
#' @param data_table data.table: data.table created from \code{create_data}
#' @param method string: selection for the spatial CRC model.
#' @param suppression int: suppression cap, defaults to 10.
#'
#' @return data.table: output of the \code{analyze} function
#' @keywords internal
analyze_wrapper <- function(x, data_table, method, suppression = 10) {
    return(
        analyze(
            data_table,
            suppress = suppression,
            method = method,
            formula_selection = "stepwise",
            opts.stepwise = c(x, verbose = FALSE)
        )
    )
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
format_analysis <- function(analyzed_dt) {
    bind_data <- rbindlist(analyzed_dt, idcol = c("Model", names(analyzed_dt)))
    return(
        bind_data[,
            Model := paste0(gsub(
                "b",
                "Backward-",
                gsub(
                    "f",
                    "Forward-",
                    gsub("fb", "Both-", Model, fixed = TRUE),
                    fixed = TRUE
                ),
                fixed = TRUE
            ))
        ]
    )
}

#' Function to run both the Poisson and Negative Binomial CRC models on the
#' same data.table and return a combined data.table.
#'
#' @param n_individuals int: number of individuals to simulate in the
#' population.
#' @param p_captures list: list of capture probabilities to emulate
#' @param p_strata list: list of extra 'stratification' probabilities to
#' emulate.
#' @param config list: list of stepwise selection parameters to run the models.
#' Each element of the list should be a list containing 'direction' and
#' 'threshold'
#'
#' @return data.table: combined data.table of results from both models
#'
#' @keywords internal
run_poisson_and_negbin <- function(
    n_individuals,
    p_captures,
    p_strata,
    config
) {
    # actually does sampling, must be done for each boostrap
    data_table <- create_data(n_individuals, p_captures, p_strata)

    pois_data <- lapply(
        config,
        analyze_wrapper,
        data_table = data_table,
        method = "poisson"
    ) |>
        format_analysis()

    nb_data <- lapply(
        config,
        analyze_wrapper,
        data_table = data_table,
        method = "negbin"
    ) |>
        format_analysis()

    return(rbind(
        nb_data[, Method := "NB"],
        pois_data[, Method := "Poisson"]
    ))
}
