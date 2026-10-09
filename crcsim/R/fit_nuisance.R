################################################################################
# File: fit_nuisance.R                                                         #
# Project: crcsim                                                              #
# Created Date: 2026-08-27                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-09-08                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

#' Estimate nuisance parameters using the specified function.
#'
#' @param func The function to use for nuisance estimation.
#' @param train The training dataset.
#' @param test The testing dataset.
#' @param n_lists The number of lists.
#' @param j Index for the first variable.
#' @param k Index for the second variable.
#' @param margin The margin for estimation.
#' @param diagnostics Controls GLM fit diagnostics.
#'
#' @return Estimated nuisance parameters.
#' @keywords internal
nuisance_estimation <- function(
    func,
    train,
    test,
    n_lists,
    j,
    k,
    margin,
    diagnostics = "quiet"
) {
    diagnostics <- match.arg(diagnostics, c("quiet", "warning", "verbose"))
    # Temporary implementation only using logit estimation
    return(qhat_logit(train, test, n_lists, j, k, margin, diagnostics))
}


#' Validate inputs shared by the nuisance estimators.
#'
#' @param train The training dataset.
#' @param test The testing dataset.
#' @param n_lists The number of lists.
#' @param j Index for the first variable.
#' @param k Index for the second variable.
#' @param margin The margin for estimation.
#' @param smooth_terms Character vector naming covariates to model with smooth
#' terms.
#' @param diagnostics The diagnostics mode.
#' @param strict Whether to apply GAM-specific data and smooth-term checks.
#' @return Validated covariate names and smooth terms.
#' @keywords internal
validate_qhat_inputs <- function(
    train,
    test,
    n_lists,
    j,
    k,
    margin,
    smooth_terms = NULL,
    diagnostics = "quiet",
    strict = FALSE
) {
    diagnostics <- match.arg(diagnostics, c("quiet", "warning", "verbose"))

    if (!is.data.frame(train) || !is.data.frame(test)) {
        stop("train and test must both be data frames.")
    }
    if (nrow(train) < 1 || nrow(test) < 1) {
        stop("train and test must both contain at least one row.")
    }
    if (
        length(n_lists) != 1 ||
            !is.numeric(n_lists) ||
            !is.finite(n_lists) ||
            n_lists < 2 ||
            n_lists %% 1 != 0
    ) {
        stop("n_lists must be an integer of at least 2.")
    }
    if (ncol(train) <= n_lists) {
        stop("train must contain at least one covariate.")
    }
    if (!identical(names(train), names(test))) {
        stop("train and test must have identical column names and order.")
    }
    if (
        length(c(j, k)) != 2 ||
            !isTRUE(is.numeric(c(j, k))) ||
            !isTRUE(all(is.finite(c(j, k)))) ||
            any(c(j, k) %% 1 != 0) ||
            any(c(j, k) < 1 | c(j, k) > n_lists) ||
            j == k
    ) {
        stop("j and k must be distinct capture-column indices.")
    }
    if (
        length(margin) != 1 ||
            !is.numeric(margin) ||
            !is.finite(margin) ||
            margin <= 0 ||
            margin >= 1
    ) {
        stop("margin must be a finite numeric value between 0 and 1.")
    }

    covariate_names <- names(train)[seq.int(n_lists + 1, ncol(train))]
    if (strict && !validate_binary_cols(train, n_lists)) {
        stop("The first n_lists columns of train must be binary.")
    }
    if (strict && !validate_binary_cols(test, n_lists)) {
        stop("The first n_lists columns of test must be binary.")
    }
    if (is.null(smooth_terms)) {
        smooth_terms <- covariate_names
    } else if (
        !is.character(smooth_terms) ||
            anyNA(smooth_terms) ||
            anyDuplicated(smooth_terms) > 0 ||
            !all(smooth_terms %in% covariate_names)
    ) {
        stop("smooth_terms must name unique covariates in train.")
    }
    if (
        strict &&
            (!isTRUE(all(vapply(
                train[covariate_names],
                is.numeric,
                logical(1L)
            ))) ||
                !isTRUE(all(vapply(
                    test[covariate_names],
                    is.numeric,
                    logical(1L)
                ))))
    ) {
        stop("All covariates must be numeric.")
    }

    return(list(
        covariate_names = covariate_names,
        smooth_terms = smooth_terms
    ))
}

#' Compress repeated nuisance-model rows and create frequency weights.
#'
#' @param fit_data The selected and consistently named modeling data.
#' @param n_rows The number of rows before compression.
#' @return A list containing compressed fit data and frequency weights.
#' @keywords internal
compress_qhat_data <- function(fit_data, n_rows) {
    if (n_rows > 10000) {
        # Expanded frequency data can contain many identical rows. Compressing
        # those rows and passing frequencies as weights preserves the fit while
        # avoiding an unnecessarily large model call.
        fit_groups <- interaction(fit_data, drop = TRUE, lex.order = TRUE)
        fit_rows <- !duplicated(fit_groups)
        fit_weights <- tabulate(
            as.integer(fit_groups),
            nbins = nlevels(fit_groups)
        )
        fit_data <- fit_data[fit_rows, , drop = FALSE]
        fit_weights <- fit_weights[as.integer(fit_groups[fit_rows])]
    } else {
        fit_weights <- rep(1, nrow(fit_data))
    }

    return(list(fit_data = fit_data, fit_weights = fit_weights))
}

#' Fit a nuisance model and record diagnostics.
#'
#' @param fit_call Function that fits and returns a nuisance model.
#' @param model_name Label used in diagnostic messages.
#' @param diagnostic_state Environment receiving diagnostic messages.
#' @return The fitted model, or NULL when fitting is unusable.
#' @keywords internal
fit_qhat_model <- function(fit_call, model_name, diagnostic_state) {
    condition_state <- new.env(parent = emptyenv())
    condition_state$warnings <- character()
    condition_state$error <- NULL
    fit <- tryCatch(
        withCallingHandlers(
            fit_call(),
            warning = function(warning) {
                condition_state$warnings <- c(
                    condition_state$warnings,
                    conditionMessage(warning)
                )
                invokeRestart("muffleWarning")
            }
        ),
        error = function(error) {
            condition_state$error <- conditionMessage(error)
            return(NULL)
        }
    )

    if (length(condition_state$warnings) > 0) {
        diagnostic_state$messages <- c(
            diagnostic_state$messages,
            paste0(model_name, ": ", unique(condition_state$warnings))
        )
    }
    if (!is.null(condition_state$error)) {
        diagnostic_state$messages <- c(
            diagnostic_state$messages,
            paste0(model_name, ": error: ", condition_state$error)
        )
    }
    if (is.null(fit)) {
        return(NULL)
    }
    if (!isTRUE(fit$converged)) {
        diagnostic_state$messages <- c(
            diagnostic_state$messages,
            paste0(model_name, ": did not converge")
        )
        return(NULL)
    }
    if (!all(is.finite(stats::coef(fit)))) {
        diagnostic_state$messages <- c(
            diagnostic_state$messages,
            paste0(model_name, ": non-finite coefficients")
        )
        return(NULL)
    }
    return(fit)
}

#' Fit a binomial GLM nuisance model.
#'
#' @param form Model formula.
#' @param model_name Label used in diagnostic messages.
#' @param fit_data Modeling data.
#' @param fit_weights Frequency weights.
#' @param max_iterations Maximum fitting iterations.
#' @param diagnostic_state Environment receiving diagnostic messages.
#' @return The fitted model, or NULL when fitting is unusable.
#' @keywords internal
fit_qhat_glm <- function(
    form,
    model_name,
    fit_data,
    fit_weights,
    max_iterations,
    diagnostic_state
) {
    fit_call <- function() {
        fit <- stats::glm(
            form,
            family = stats::binomial(link = "logit"),
            data = fit_data,
            weights = fit_weights,
            control = stats::glm.control(maxit = max_iterations)
        )
        return(fit)
    }
    fit <- fit_qhat_model(
        fit_call,
        model_name,
        diagnostic_state
    )
    return(fit)
}

#' Fit a binomial GAM nuisance model.
#'
#' @param form Model formula.
#' @param model_name Label used in diagnostic messages.
#' @param fit_data Modeling data.
#' @param fit_weights Frequency weights.
#' @param diagnostic_state Environment receiving diagnostic messages.
#' @return The fitted model, or NULL when fitting is unusable.
#' @keywords internal
fit_qhat_gam <- function(
    form,
    model_name,
    fit_data,
    fit_weights,
    diagnostic_state
) {
    fit_call <- function() {
        fit <- mgcv::gam(
            form,
            family = stats::binomial(link = "logit"),
            data = fit_data,
            weights = fit_weights
        )
        return(fit)
    }
    fit <- fit_qhat_model(
        fit_call,
        model_name,
        diagnostic_state
    )
    return(fit)
}

#' Predict nuisance probabilities and record diagnostics.
#'
#' @param predict_call Function that returns model predictions.
#' @param model_name Label used in diagnostic messages.
#' @param n_test Number of expected predictions.
#' @param margin Minimum allowed probability.
#' @param diagnostic_state Environment receiving diagnostic messages.
#' @return Bounded predictions, or NULL when prediction is unusable.
#' @keywords internal
predict_qhat_model <- function(
    predict_call,
    model_name,
    n_test,
    margin,
    diagnostic_state
) {
    condition_state <- new.env(parent = emptyenv())
    condition_state$warnings <- character()
    condition_state$error <- NULL
    predictions <- tryCatch(
        withCallingHandlers(
            predict_call(),
            warning = function(warning) {
                condition_state$warnings <- c(
                    condition_state$warnings,
                    conditionMessage(warning)
                )
                invokeRestart("muffleWarning")
            }
        ),
        error = function(error) {
            condition_state$error <- conditionMessage(error)
            return(NULL)
        }
    )

    if (length(condition_state$warnings) > 0) {
        diagnostic_state$messages <- c(
            diagnostic_state$messages,
            paste0(model_name, ": ", unique(condition_state$warnings))
        )
    }
    if (!is.null(condition_state$error)) {
        diagnostic_state$messages <- c(
            diagnostic_state$messages,
            paste0(model_name, ": error: ", condition_state$error)
        )
    }
    if (is.null(predictions)) {
        return(NULL)
    }
    if (length(predictions) != n_test) {
        diagnostic_state$messages <- c(
            diagnostic_state$messages,
            paste0(model_name, ": unexpected prediction length")
        )
        return(NULL)
    }
    if (!isTRUE(all(is.finite(predictions)))) {
        diagnostic_state$messages <- c(
            diagnostic_state$messages,
            paste0(model_name, ": non-finite predictions")
        )
        return(NULL)
    }
    return(pmin(pmax(as.numeric(predictions), margin), 1))
}

#' Predict from a GAM nuisance model.
#'
#' @param fit Fitted GAM model.
#' @param test Testing data.
#' @param model_name Label used in diagnostic messages.
#' @param margin Minimum allowed probability.
#' @param diagnostic_state Environment receiving diagnostic messages.
#' @return Bounded predictions, or NULL when prediction is unusable.
#' @keywords internal
predict_qhat_gam <- function(
    fit,
    test,
    model_name,
    margin,
    diagnostic_state
) {
    predict_call <- function() {
        predictions <- stats::predict(fit, newdata = test, type = "response")
        return(predictions)
    }
    predictions <- predict_qhat_model(
        predict_call,
        model_name,
        nrow(test),
        margin,
        diagnostic_state
    )
    return(predictions)
}

#' Estimate the initial nuisance parameters using GLM fits on the training data
#' and predict on the test data.
#'
#' @param train The training dataset.
#' @param test The testing dataset.
#' @param n_lists The number of lists.
#' @param j Index for the first variable.
#' @param k Index for the second variable.
#' @param margin The margin for estimation.
#' @param diagnostics Controls GLM fit diagnostics.
#'
#' @importFrom stats predict glm binomial
#' @return A list containing the initial estimates for q_1, q_2, and q_12.
#' @keywords internal
qhat_logit <- function(
    train,
    test,
    n_lists,
    j,
    k,
    margin,
    diagnostics = "quiet"
) {
    validate_qhat_inputs(
        train,
        test,
        n_lists,
        j,
        k,
        margin,
        diagnostics = diagnostics
    )
    model_columns <- c(j, k, (n_lists + 1):ncol(train))
    fit_data <- train[, model_columns, drop = FALSE]
    model_names <- c(
        paste0("d", seq_len(n_lists)),
        paste0("x", seq_len(ncol(train) - n_lists))
    )
    colnames(fit_data) <- model_names[c(j, k, (n_lists + 1):ncol(train))]

    compressed_data <- compress_qhat_data(fit_data, nrow(train))
    fit_data <- compressed_data$fit_data
    fit_weights <- compressed_data$fit_weights

    diagnostic_state <- new.env(parent = emptyenv())
    diagnostic_state$messages <- character()

    template_q_j <- function(d_fit, q_12_offset = 0.0) {
        return(
            pmax(
                q_12_offset + predict(d_fit, newdata = test, type = "response"),
                margin
            )
        )
    }

    ## Core Functionality Start
    c_names <- c(paste0("d", 1:n_lists), paste0("x", 1:(ncol(train) - n_lists)))
    colnames(train) <- c_names
    colnames(test) <- c_names

    max_iterations <- if (length(fit_weights) < nrow(train)) 1000 else 25
    fit_j_0 <- fit_qhat_glm(
        formula(paste0("d", j, "*(1 - d", k, ") ~.")),
        paste0("d", j, "*(1-d", k, ")"),
        fit_data,
        fit_weights,
        max_iterations,
        diagnostic_state
    )
    fit_0_k <- fit_qhat_glm(
        formula(paste0("d", k, "*(1 - d", j, ") ~.")),
        paste0("d", k, "*(1-d", j, ")"),
        fit_data,
        fit_weights,
        max_iterations,
        diagnostic_state
    )
    fit_j_k <- fit_qhat_glm(
        formula(paste0("d", j, "*d", k, " ~.")),
        paste0("d", j, "*d", k),
        fit_data,
        fit_weights,
        max_iterations,
        diagnostic_state
    )

    if (length(diagnostic_state$messages) > 0) {
        diagnostic_message <- paste(
            c(
                "GLM nuisance diagnostics:",
                paste0("- ", diagnostic_state$messages)
            ),
            collapse = "\n"
        )
        if (diagnostics == "warning") {
            warning(diagnostic_message, call. = FALSE)
        } else if (diagnostics == "verbose") {
            message(diagnostic_message)
        }
    }

    if (
        is.null(fit_j_0) ||
            is.null(fit_0_k) ||
            is.null(fit_j_k)
    ) {
        return(NULL)
    }

    q_12 <- template_q_j(fit_j_k)
    q_1 <- pmin(pmax(q_12, template_q_j(fit_j_0, q_12_offset = q_12)), 1)
    q_2 <- pmax(
        q_12 / q_1,
        pmin(template_q_j(fit_0_k, q_12_offset = q_12), 1 + q_12 - q_1, 1)
    )

    return(list(q_1 = q_1, q_2 = q_2, q_12 = q_12))
}

#' Estimate initial nuisance parameters using GAM fits.
#'
#' This is an implementation outline only. The current function deliberately
#' remains unavailable until the smooth-term interface is agreed upon.
#'
#' @param train The training dataset.
#' @param test The testing dataset.
#' @param n_lists The number of lists.
#' @param j Index for the first variable.
#' @param k Index for the second variable.
#' @param margin The margin for estimation.
#' @param smooth_terms Character vector naming covariates to model with smooth
#' terms. If NULL, all eligible covariates are smoothed by default.
#' @param diagnostics Controls GAM fit diagnostics.
#' @return A list containing the initial estimates for q_1, q_2, and q_12.
#' @keywords internal
qhat_gam <- function(
    train,
    test,
    n_lists,
    j,
    k,
    margin,
    smooth_terms = NULL,
    diagnostics = "quiet"
) {
    input_info <- validate_qhat_inputs(
        train,
        test,
        n_lists,
        j,
        k,
        margin,
        smooth_terms,
        diagnostics,
        strict = TRUE
    )
    smooth_terms <- input_info$smooth_terms

    model_columns <- c(j, k, (n_lists + 1):ncol(train))
    model_names <- c(
        paste0("d", seq_len(n_lists)),
        paste0("x", seq_len(ncol(train) - n_lists))
    )
    fit_data <- train[, model_columns, drop = FALSE]
    colnames(fit_data) <- model_names[model_columns]

    colnames(train) <- model_names
    colnames(test) <- model_names
    smooth_terms <- paste0(
        "x",
        match(smooth_terms, input_info$covariate_names)
    )

    compressed_data <- compress_qhat_data(fit_data, nrow(train))
    fit_data <- compressed_data$fit_data
    fit_weights <- compressed_data$fit_weights

    response_j <- paste0("d", j)
    response_k <- paste0("d", k)
    fit_data$q_j0 <- fit_data[[response_j]] * (1 - fit_data[[response_k]])
    fit_data$q_0k <- (1 - fit_data[[response_j]]) * fit_data[[response_k]]
    fit_data$q_jk <- fit_data[[response_j]] * fit_data[[response_k]]

    covariate_names <- model_names[(n_lists + 1):ncol(train)]
    linear_terms <- setdiff(covariate_names, smooth_terms)
    smooth_formula_terms <- paste0("s(", smooth_terms, ")")
    formula_terms <- c(smooth_formula_terms, linear_terms)
    formula_rhs <- if (length(formula_terms) == 0) {
        "1"
    } else {
        paste(formula_terms, collapse = " + ")
    }
    gam_formulas <- list(
        q_j0 = stats::as.formula(paste("q_j0 ~", formula_rhs)),
        q_0k = stats::as.formula(paste("q_0k ~", formula_rhs)),
        q_jk = stats::as.formula(paste("q_jk ~", formula_rhs))
    )

    diagnostic_state <- new.env(parent = emptyenv())
    diagnostic_state$messages <- character()

    fit_j_0 <- fit_qhat_gam(
        gam_formulas$q_j0,
        "q_j0",
        fit_data,
        fit_weights,
        diagnostic_state
    )
    fit_0_k <- fit_qhat_gam(
        gam_formulas$q_0k,
        "q_0k",
        fit_data,
        fit_weights,
        diagnostic_state
    )
    fit_j_k <- fit_qhat_gam(
        gam_formulas$q_jk,
        "q_jk",
        fit_data,
        fit_weights,
        diagnostic_state
    )

    if (length(diagnostic_state$messages) > 0) {
        diagnostic_message <- paste(
            c(
                "GAM nuisance diagnostics:",
                paste0("- ", diagnostic_state$messages)
            ),
            collapse = "\n"
        )
        if (diagnostics == "warning") {
            warning(diagnostic_message, call. = FALSE)
        } else if (diagnostics == "verbose") {
            message(diagnostic_message)
        }
    }

    if (
        is.null(fit_j_0) ||
            is.null(fit_0_k) ||
            is.null(fit_j_k)
    ) {
        return(NULL)
    }

    prediction_message_start <- length(diagnostic_state$messages)
    q_j0_hat <- predict_qhat_gam(
        fit_j_0,
        test,
        "q_j0",
        margin,
        diagnostic_state
    )
    q_0k_hat <- predict_qhat_gam(
        fit_0_k,
        test,
        "q_0k",
        margin,
        diagnostic_state
    )
    q_jk_hat <- predict_qhat_gam(
        fit_j_k,
        test,
        "q_jk",
        margin,
        diagnostic_state
    )

    prediction_messages <- diagnostic_state$messages[
        seq.int(prediction_message_start + 1, length(diagnostic_state$messages))
    ]
    if (length(prediction_messages) > 0) {
        diagnostic_message <- paste(
            c(
                "GAM nuisance diagnostics:",
                paste0("- ", prediction_messages)
            ),
            collapse = "\n"
        )
        if (diagnostics == "warning") {
            warning(diagnostic_message, call. = FALSE)
        } else if (diagnostics == "verbose") {
            message(diagnostic_message)
        }
    }

    if (
        is.null(q_j0_hat) ||
            is.null(q_0k_hat) ||
            is.null(q_jk_hat)
    ) {
        return(NULL)
    }

    # 7. Combine the predicted exclusive-pattern probabilities with the
    #    intersection and enforce the same compatibility constraints as the
    #    GLM nuisance estimates.
    q_12 <- q_jk_hat
    q_1 <- pmin(pmax(q_12, q_12 + q_j0_hat), 1)
    q_2 <- pmax(
        q_12 / q_1,
        pmin(q_12 + q_0k_hat, 1 + q_12 - q_1, 1)
    )

    # 8. Return the bounded initial nuisance estimates. Unusable fits or
    #    predictions return NULL above, after diagnostics are reported.
    return(list(q_1 = q_1, q_2 = q_2, q_12 = q_12))
}

#' Fit TMLE nuisance parameters. This updates the initial estimates using the
#' TMLE procedure.
#'
#' @param q_1 Initial estimate for the first nuisance parameter.
#' @param q_2 Initial estimate for the second nuisance parameter.
#' @param q_12 Initial estimate for the joint nuisance parameter.
#' @param y_j Observed outcome for the first variable.
#' @param y_k Observed outcome for the second variable.
#' @param y_jk Observed joint outcome for the first and second variables.
#' @param iterations Maximum number of iterations for the TMLE update.
#' @param margin Margin for estimation.
#' @param n_lists Number of lists.
#'
#' @return Updated nuisance parameter estimates.
#' @keywords internal
tmle_nuisance <- function(
    q_1,
    q_2,
    q_12,
    y_j,
    y_k,
    y_jk,
    iterations,
    margin,
    n_lists
) {
    ############################################################################
    # Helper Functions
    ############################################################################
    logit <- function(x) {
        return(log(x / (1 - x)))
    }

    ############################################################################
    # Reset Nuisance Values
    ############################################################################

    margin_error <- 1 + margin
    count <- 0

    q_10 <- pmin(pmax(q_1 - q_12, margin), 1 - margin)
    q_02 <- pmin(pmax(q_2 - q_12, margin), 1 - margin)
    q_12 <- pmin(pmax(q_12, margin), 1 - margin)

    y_j0 <- y_j * (1 - y_k)
    y_0k <- (1 - y_j) * y_k

    while (abs(margin_error) > margin && count < iterations) {
        q12_fit <- fit_glm_tmle(
            y_jk,
            logit(q_12),
            (q_10 + q_12) /
                q_12 +
                (q_02 + q_12) / q_12 -
                (q_10 + q_12) * (q_02 + q_12) / q_12^2
        )
        if (!is.null(q12_fit$value)) {
            q_12 <- pmax(pmin(q12_fit$value, 1), margin)
        }

        q10_fit <- fit_glm_tmle(
            y_j0,
            logit(q_10),
            (q_02 + q_12) / q_12
        )
        if (!is.null(q10_fit$value)) {
            q_10 <- pmax(pmin(q10_fit$value, 1 - q_12), margin)
        }

        if (n_lists > 2) {
            q02_fit <- fit_glm_tmle(
                y_0k,
                logit(q_02),
                (q_10 + q_12) / q_12
            )
            if (!is.null(q02_fit$value)) {
                q_02 <- pmax(
                    pmin(q02_fit$value, 1 - q_10 - q_12),
                    margin
                )
            }
        } else {
            q02_fit <- list(error = 0)
            q_02 <- pmax(0, 1 - q_10 - q_12)
        }

        margin_error <- max(q12_fit$error, q10_fit$error, q02_fit$error)
        count <- count + 1
    }

    new_nuisances <- list(
        q_1 = q_10 + q_12,
        q_2 = q_02 + q_12,
        q_12 = q_12
    )
    return(new_nuisances)
}

#' Fit a GLM for TMLE updates.
#'
#' @param response The response variable.
#' @param offset The offset for the GLM.
#' @param ratio The ratio used as the predictor in the GLM.
#'
#' @return A list containing the fitted values and the maximum absolute
#' coefficient value as the error.
#'
#' @importFrom stats complete.cases glm.fit coef binomial
#' @keywords internal
fit_glm_tmle <- function(response, offset, ratio) {
    complete <- complete.cases(response, offset, ratio)

    fit <- try(
        glm.fit(
            x = matrix(ratio[complete], ncol = 1),
            y = response[complete],
            offset = offset[complete],
            family = binomial(link = "logit")
        ),
        silent = TRUE
    )

    if (inherits(fit, "try-error")) {
        return(list(value = NULL, error = Inf))
    }

    return(list(
        value = fit$fitted.values,
        error = max(abs(coef(fit)), na.rm = TRUE)
    ))
}
