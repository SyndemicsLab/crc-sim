################################################################################
# File: covariate_control.R                                                    #
# Project: crcsim                                                              #
# Created Date: 2026-07-29                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-08-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Build Covariate Data
#' @description This function generates a data frame containing binary
#' covariate data for a specified number of individuals and covariates. Each
#' covariate is generated as a binary variable (0 or 1) with a probability of
#' 0.5.
#'
#' @param n_individuals The number of individuals for whom to generate
#' covariate data.
#' @param n_covariates The number of covariates to generate for each individual.
#' @return A data frame with n_individuals rows and n_covariates columns, where
#' each column represents a covariate and each row represents an individual.
#'
#' importFrom stats rbinom setNames
#' @export
build_covariate_data <- function(n_individuals, n_covariates) {
    if (n_covariates == 0) {
        return(data.frame(row.names = seq_len(n_individuals)))
    }

    cov_data <- data.frame(
        matrix(
            rbinom(n_individuals * n_covariates, 1, 0.5),
            nrow = n_individuals,
            ncol = n_covariates
        )
    ) |>
        setNames(paste0("covariate_", seq_len(n_covariates)))
    return(cov_data)
}

#' Relocate Covariates
#' @description This function relocates specified covariate columns in a data
#' frame to the end of the data frame. It takes a data frame and a vector of
#' covariate column names as input and returns the modified data frame with the
#' covariate columns relocated to the end.
#'
#' @param data A data frame containing the covariate columns to be relocated
#' @param covariate_cols A character vector of covariate column names to be
#' relocated.
#' @return A data frame with the specified covariate columns relocated to the
#' end.
#'
#' @importFrom dplyr relocate last_col
#' @importFrom rlang sym
#' @export
relocate_covariates <- function(data, covariate_cols) {
    for (col in covariate_cols) {
        data <- data |>
            relocate(!!sym(col), .after = last_col())
    }
    return(data)
}
