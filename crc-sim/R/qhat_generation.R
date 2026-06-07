################################################################################
# File: qhat_generation.R                                                      #
# Project: crc-sim                                                             #
# Created Date: 2026-05-15                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-28                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Generate nuisance function estimates for one train/test split
#'
#' Fits a single nuisance function on a training split and predicts capture
#' probabilities on a test split for all pairwise capture-list combinations.
#'
#' @param train Training data frame with capture columns followed by covariates.
#' @param test Test data frame with the same columns as train.
#' @param capture_indicators Character vector of capture column names, or
#'   integer vector of column indices (1-indexed).
#' @param nuisance_function One of "logit", "gam", or "ranger".
#' @param margin Numeric threshold used for overlap warnings and nuisance
#'   truncation. Default is 0.005.
#' @param ... Additional arguments passed to the qhat_* function.
#'
#' @return A nested list keyed by pair names. Each element contains a list with
#'   q1, q2, and q12 predictions for the rows in test.
#'
#' @keywords internal
qhat_generation <- function(
    train,
    test,
    j,
    k,
    capture_names,
    nuisance_function = c("glm", "gam", "random_forest"),
    margin = 0.005,
    ...
) {
    nuisance_function <- match.arg(nuisance_function)

    if (!is.data.frame(train) || !is.data.frame(test)) {
        stop("`train` and `test` must both be data frames.")
    }

    if (!identical(names(train), names(test))) {
        stop("`train` and `test` must have identical column names.")
    }

    capture_idx <- match(capture_names, names(train))
    if (anyNA(capture_idx)) {
        stop("Some `capture_names` were not found in the data.")
    }

    invalid_train <- !vapply(train[capture_idx], is_column_binary, logical(1))
    invalid_test <- !vapply(test[capture_idx], is_column_binary, logical(1))
    if (any(invalid_train) || any(invalid_test)) {
        invalid_names <- unique(c(
            capture_names[invalid_train],
            capture_names[invalid_test]
        ))
        stop(
            paste0(
                "Capture columns must be binary in both train and test: ",
                toString(invalid_names),
                "."
            )
        )
    }

    estimate_qhats <- get(
        paste0(nuisance_function, "_qhats"),
        mode = "function"
    )

    overlap <- mean(train[[j]] * train[[k]], na.rm = TRUE)
    if (overlap < margin) {
        warning(
            paste0(
                "Overlap between capture columns %s and %s is ",
                "%.6f, below margin %.6f."
            ),
            j,
            k,
            overlap,
            margin
        )
    }

    qhat_raw <- estimate_qhats(
        train = train,
        test = test,
        n_capture_cols = length(capture_idx),
        j = j,
        k = k,
        margin = margin,
        ...
    )

    # Guarantee q1 is less than or equal to 1 and q2 is greater than or
    # equal to q12 / q1, using the fact that q12 <= q1 and q12 <= q2. This
    # is a post-processing step to ensure valid probabilities.
    q12 <- qhat_raw[[q12]]
    q1 <- pmin(pmax(q12, qhat_raw[[q1]]), 1)
    q2 <- pmax(q12 / q1, pmin(qhat_raw[[q2]], 1 + q12 - q1, 1))

    return(list(
        q1 = q1,
        q2 = q2,
        q12 = q12,
        y_j = test[[j]],
        y_k = test[[k]]
    ))
}
