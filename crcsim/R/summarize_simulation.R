#' Summarize simulated or observed capture-recapture data
#'
#' Create reusable capture-history, list-overlap, and covariate summaries from
#' a structured simulation or a raw data frame. Raw data require an explicit
#' scope because the all-zero population is not identifiable from observed-only
#' data.
#'
#' @param data A \\code{crcsim_simulation} object or a data frame containing
#'   binary capture columns.
#' @param capture_columns Character vector naming the capture columns.
#' @param covariate_columns Optional character vector naming covariate columns.
#'   When omitted, non-capture columns are treated as covariates.
#' @param scope Required for raw data. One of \\code{full},
#'   \\code{observed}, or \\code{unobserved}.
#' @details
#' A structured simulation supplies all three scopes and its known generating
#' truth. Raw data supply only the explicitly declared scope, so the other
#' scope rows are retained with missing values rather than inferred. The
#' \\code{capture_probabilities} table summarizes known full-population truth;
#' it is unavailable for raw data.
#' @returns A classed list containing tidy diagnostic tables and metadata.
#'   Numeric covariates are summarized in \\code{covariates}; categorical
#'   covariates are summarized in \\code{covariate_levels}; known true
#'   probabilities are summarized in \\code{capture_probabilities}.
#'
#' @export
summarize_simulation <- function(
    data,
    capture_columns,
    covariate_columns = NULL,
    scope = NULL
) {
    canonical_scopes <- c(
        "full",
        "observed",
        "unobserved"
    )
    input <- normalize_summary_input(
        data = data,
        scope = scope,
        canonical_scopes = canonical_scopes
    )
    validate_summary_columns(
        views = input$views,
        capture_columns = capture_columns
    )

    if (is.null(covariate_columns)) {
        covariate_columns <- setdiff(
            names(input$views[[input$available_scopes[[1]]]]),
            capture_columns
        )
    }
    validate_covariate_columns(input$views, covariate_columns)

    counts <- summarize_scope_counts(
        views = input$views,
        available_scopes = input$available_scopes,
        canonical_scopes = canonical_scopes
    )
    capture_patterns <- summarize_capture_patterns(
        views = input$views,
        capture_columns = capture_columns,
        available_scopes = input$available_scopes,
        canonical_scopes = canonical_scopes
    )
    list_counts <- summarize_list_counts(
        views = input$views,
        capture_columns = capture_columns,
        available_scopes = input$available_scopes,
        canonical_scopes = canonical_scopes
    )
    overlaps <- summarize_overlaps(
        views = input$views,
        capture_columns = capture_columns,
        available_scopes = input$available_scopes,
        canonical_scopes = canonical_scopes
    )
    covariates <- summarize_covariates(
        views = input$views,
        covariate_columns = covariate_columns,
        available_scopes = input$available_scopes,
        canonical_scopes = canonical_scopes
    )
    covariate_levels <- summarize_covariate_levels(
        views = input$views,
        covariate_columns = covariate_columns,
        available_scopes = input$available_scopes,
        canonical_scopes = canonical_scopes
    )
    capture_probabilities <- summarize_capture_probabilities(
        truth = input$truth,
        capture_columns = capture_columns,
        canonical_scopes = canonical_scopes
    )

    metadata <- input$metadata
    metadata$summary_version <- "0.1.0"
    metadata$input_type <- input$input_type
    metadata$available_scopes <- input$available_scopes
    metadata$capture_columns <- capture_columns
    metadata$covariate_columns <- covariate_columns

    result <- list(
        counts = counts,
        capture_patterns = capture_patterns,
        list_counts = list_counts,
        overlaps = overlaps,
        covariates = covariates,
        covariate_levels = covariate_levels,
        capture_probabilities = capture_probabilities,
        metadata = metadata
    )
    class(result) <- c("crcsim_summary", "list")
    return(result)
}

normalize_summary_input <- function(data, scope, canonical_scopes) {
    if (inherits(data, "crcsim_simulation")) {
        views <- list(
            full = data$full_data,
            observed = data$observed_data,
            unobserved = data$unobserved_data
        )
        available_scopes <- canonical_scopes
        metadata <- data$metadata
        truth <- data$truth
        input_type <- "crcsim_simulation"
    } else if (is.data.frame(data)) {
        if (
            is.null(scope) || length(scope) != 1 || !scope %in% canonical_scopes
        ) {
            stop("scope is required for raw data and must be a valid scope")
        }
        views <- setNames(list(tibble::as_tibble(data)), scope)
        available_scopes <- scope
        metadata <- list()
        truth <- NULL
        input_type <- "data.frame"
    } else {
        stop("data must be a crcsim_simulation or data.frame")
    }

    return(list(
        views = views,
        available_scopes = available_scopes,
        metadata = metadata,
        truth = truth,
        input_type = input_type
    ))
}

validate_summary_columns <- function(views, capture_columns) {
    if (!is.character(capture_columns) || length(capture_columns) < 1) {
        stop("capture_columns must be a non-empty character vector")
    }
    for (view in views) {
        if (!all(capture_columns %in% names(view))) {
            stop("capture_columns must name columns in every input view")
        }
        if (
            nrow(view) > 0 &&
                !validate_binary_cols(
                    view,
                    end = match(
                        capture_columns[length(capture_columns)],
                        names(view)
                    ),
                    start = match(capture_columns[1], names(view))
                )
        ) {
            if (
                any(
                    vapply(
                        view[capture_columns],
                        function(column) {
                            all(column %in% c(0, 1))
                        },
                        logical(1)
                    ) ==
                        FALSE
                )
            ) {
                stop("capture columns must contain only 0 and 1")
            }
        }
    }
    return(invisible(NULL))
}

validate_covariate_columns <- function(views, covariate_columns) {
    if (!is.character(covariate_columns)) {
        stop("covariate_columns must be a character vector or NULL")
    }
    for (view in views) {
        if (!all(covariate_columns %in% names(view))) {
            stop("covariate_columns must name columns in every input view")
        }
    }
    return(invisible(NULL))
}

add_unavailable_scopes <- function(summary, canonical_scopes, columns) {
    missing_scopes <- setdiff(canonical_scopes, unique(summary$scope))
    if (length(missing_scopes) == 0) {
        return(summary)
    }
    unavailable <- tibble::tibble(scope = missing_scopes)
    for (column in columns) {
        unavailable[[column]] <- NA
    }
    return(dplyr::bind_rows(summary, unavailable))
}

summarize_scope_counts <- function(views, available_scopes, canonical_scopes) {
    result <- dplyr::bind_rows(lapply(available_scopes, function(scope) {
        tibble::tibble(
            scope = scope,
            n_individuals = nrow(views[[scope]])
        )
    }))
    full_count <- result$n_individuals[result$scope == "full"]
    result$proportion_of_full <- if (length(full_count) == 1) {
        result$n_individuals / full_count
    } else {
        NA_real_
    }
    return(add_unavailable_scopes(
        result,
        canonical_scopes,
        c("n_individuals", "proportion_of_full")
    ))
}

summarize_capture_patterns <- function(
    views,
    capture_columns,
    available_scopes,
    canonical_scopes
) {
    result <- dplyr::bind_rows(lapply(available_scopes, function(scope) {
        view <- views[[scope]]
        if (nrow(view) == 0) {
            return(tibble::tibble(
                scope = scope,
                history = paste(
                    rep("0", length(capture_columns)),
                    collapse = ""
                ),
                n = 0,
                proportion = NA_real_
            ))
        }
        histories <- apply(view[capture_columns], 1, paste0, collapse = "")
        dplyr::count(
            tibble::tibble(history = histories),
            history,
            name = "n"
        ) |>
            dplyr::mutate(
                scope = scope,
                proportion = n / sum(n),
                .before = 1
            )
    }))
    return(add_unavailable_scopes(
        result,
        canonical_scopes,
        c("history", "n", "proportion")
    ))
}

summarize_list_counts <- function(
    views,
    capture_columns,
    available_scopes,
    canonical_scopes
) {
    result <- dplyr::bind_rows(lapply(available_scopes, function(scope) {
        view <- views[[scope]]
        dplyr::bind_rows(lapply(capture_columns, function(capture_column) {
            n <- sum(view[[capture_column]] == 1)
            tibble::tibble(
                scope = scope,
                capture = capture_column,
                n = n,
                proportion = if (nrow(view) > 0) n / nrow(view) else NA_real_
            )
        }))
    }))
    return(add_unavailable_scopes(
        result,
        canonical_scopes,
        c("capture", "n", "proportion")
    ))
}

summarize_overlaps <- function(
    views,
    capture_columns,
    available_scopes,
    canonical_scopes
) {
    if (length(capture_columns) < 2) {
        result <- dplyr::bind_rows(lapply(available_scopes, function(scope) {
            tibble::tibble(
                scope = scope,
                capture_1 = NA_character_,
                capture_2 = NA_character_,
                intersection = NA_real_,
                union = NA_real_,
                jaccard = NA_real_,
                expected_intersection = NA_real_,
                observed_to_expected = NA_real_
            )
        }))
        return(add_unavailable_scopes(
            result,
            canonical_scopes,
            c(
                "capture_1",
                "capture_2",
                "intersection",
                "union",
                "jaccard",
                "expected_intersection",
                "observed_to_expected"
            )
        ))
    }
    pairs <- utils::combn(capture_columns, 2, simplify = FALSE)
    result <- dplyr::bind_rows(lapply(available_scopes, function(scope) {
        view <- views[[scope]]
        dplyr::bind_rows(lapply(pairs, function(pair) {
            first <- view[[pair[[1]]]] == 1
            second <- view[[pair[[2]]]] == 1
            intersection <- sum(first & second)
            union <- sum(first | second)
            expected <- if (nrow(view) > 0) {
                sum(first) * sum(second) / nrow(view)
            } else {
                NA_real_
            }
            tibble::tibble(
                scope = scope,
                capture_1 = pair[[1]],
                capture_2 = pair[[2]],
                intersection = intersection,
                union = union,
                jaccard = if (union > 0) intersection / union else NA_real_,
                expected_intersection = expected,
                observed_to_expected = if (!is.na(expected) && expected > 0) {
                    intersection / expected
                } else {
                    NA_real_
                }
            )
        }))
    }))
    return(add_unavailable_scopes(
        result,
        canonical_scopes,
        c(
            "capture_1",
            "capture_2",
            "intersection",
            "union",
            "jaccard",
            "expected_intersection",
            "observed_to_expected"
        )
    ))
}

summarize_covariates <- function(
    views,
    covariate_columns,
    available_scopes,
    canonical_scopes
) {
    if (length(covariate_columns) == 0) {
        return(tibble::tibble(
            scope = character(),
            covariate = character(),
            n = integer(),
            mean = numeric(),
            sd = numeric(),
            min = numeric(),
            q25 = numeric(),
            median = numeric(),
            q75 = numeric(),
            max = numeric()
        ))
    }
    result <- dplyr::bind_rows(lapply(available_scopes, function(scope) {
        view <- views[[scope]]
        dplyr::bind_rows(lapply(covariate_columns, function(covariate) {
            values <- view[[covariate]]
            numeric_values <- if (is.numeric(values)) {
                as.numeric(values)
            } else {
                rep(NA_real_, length(values))
            }
            n_values <- sum(!is.na(numeric_values))
            tibble::tibble(
                scope = scope,
                covariate = covariate,
                n = n_values,
                mean = if (n_values > 0) {
                    mean(numeric_values, na.rm = TRUE)
                } else {
                    NA_real_
                },
                sd = if (n_values > 1) {
                    stats::sd(numeric_values, na.rm = TRUE)
                } else {
                    NA_real_
                },
                min = if (n_values > 0) {
                    min(numeric_values, na.rm = TRUE)
                } else {
                    NA_real_
                },
                q25 = if (n_values > 0) {
                    stats::quantile(
                        numeric_values,
                        0.25,
                        na.rm = TRUE,
                        names = FALSE
                    )
                } else {
                    NA_real_
                },
                median = if (n_values > 0) {
                    stats::median(numeric_values, na.rm = TRUE)
                } else {
                    NA_real_
                },
                q75 = if (n_values > 0) {
                    stats::quantile(
                        numeric_values,
                        0.75,
                        na.rm = TRUE,
                        names = FALSE
                    )
                } else {
                    NA_real_
                },
                max = if (n_values > 0) {
                    max(numeric_values, na.rm = TRUE)
                } else {
                    NA_real_
                }
            )
        }))
    }))
    return(add_unavailable_scopes(
        result,
        canonical_scopes,
        c("covariate", "n", "mean", "sd", "min", "q25", "median", "q75", "max")
    ))
}

summarize_covariate_levels <- function(
    views,
    covariate_columns,
    available_scopes,
    canonical_scopes
) {
    categorical_columns <- covariate_columns[vapply(
        views[[available_scopes[[1]]]][covariate_columns],
        function(column) !is.numeric(column),
        logical(1)
    )]
    if (length(categorical_columns) == 0) {
        return(tibble::tibble(
            scope = character(),
            covariate = character(),
            level = character(),
            n = integer(),
            proportion = numeric()
        ))
    }

    result <- dplyr::bind_rows(lapply(available_scopes, function(scope) {
        view <- views[[scope]]
        dplyr::bind_rows(lapply(categorical_columns, function(covariate) {
            levels <- unique(unlist(lapply(
                views[available_scopes],
                function(data) {
                    values <- as.character(data[[covariate]])
                    values[is.na(values)] <- "(Missing)"
                    return(values)
                }
            )))
            values <- as.character(view[[covariate]])
            values[is.na(values)] <- "(Missing)"
            dplyr::bind_rows(lapply(levels, function(level) {
                n <- sum(values == level)
                tibble::tibble(
                    scope = scope,
                    covariate = covariate,
                    level = level,
                    n = n,
                    proportion = if (length(values) > 0) {
                        n / length(values)
                    } else {
                        NA_real_
                    }
                )
            }))
        }))
    }))
    return(add_unavailable_scopes(
        result,
        canonical_scopes,
        c("covariate", "level", "n", "proportion")
    ))
}

summarize_capture_probabilities <- function(
    truth,
    capture_columns,
    canonical_scopes
) {
    probabilities <- NULL
    if (!is.null(truth)) {
        probabilities <- truth$capture_probabilities
        if (is.null(probabilities) && !is.null(truth$p_captures)) {
            probabilities <- matrix(
                rep(truth$p_captures, truth$n_total),
                nrow = truth$n_total,
                byrow = TRUE
            )
            colnames(probabilities) <- capture_columns
        }
    }
    if (is.null(probabilities)) {
        result <- tibble::tibble(
            scope = character(),
            capture = character(),
            n = integer(),
            mean = numeric(),
            sd = numeric(),
            min = numeric(),
            q25 = numeric(),
            median = numeric(),
            q75 = numeric(),
            max = numeric()
        )
        return(add_unavailable_probability_scopes(
            result,
            canonical_scopes,
            capture_columns
        ))
    }

    result <- dplyr::bind_rows(lapply(
        seq_along(capture_columns),
        function(index) {
            values <- as.numeric(probabilities[, index])
            tibble::tibble(
                scope = "full",
                capture = capture_columns[[index]],
                n = sum(!is.na(values)),
                mean = mean(values, na.rm = TRUE),
                sd = if (sum(!is.na(values)) > 1) {
                    stats::sd(values, na.rm = TRUE)
                } else {
                    NA_real_
                },
                min = min(values, na.rm = TRUE),
                q25 = stats::quantile(
                    values,
                    0.25,
                    na.rm = TRUE,
                    names = FALSE
                ),
                median = stats::median(values, na.rm = TRUE),
                q75 = stats::quantile(
                    values,
                    0.75,
                    na.rm = TRUE,
                    names = FALSE
                ),
                max = max(values, na.rm = TRUE)
            )
        }
    ))
    return(add_unavailable_probability_scopes(
        result,
        canonical_scopes,
        capture_columns
    ))
}

add_unavailable_probability_scopes <- function(
    summary,
    canonical_scopes,
    capture_columns
) {
    missing_scopes <- setdiff(canonical_scopes, unique(summary$scope))
    if (length(missing_scopes) == 0) {
        return(summary)
    }
    unavailable <- dplyr::bind_rows(lapply(missing_scopes, function(scope) {
        tibble::tibble(
            scope = scope,
            capture = capture_columns,
            n = NA_integer_,
            mean = NA_real_,
            sd = NA_real_,
            min = NA_real_,
            q25 = NA_real_,
            median = NA_real_,
            q75 = NA_real_,
            max = NA_real_
        )
    }))
    return(dplyr::bind_rows(summary, unavailable))
}
