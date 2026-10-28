test_that("estimate_no_covariates returns pairwise population estimates", {
    data <- data.frame(
        capture_1 = c(1, 0, 1, 1),
        capture_2 = c(0, 1, 1, 1)
    )

    result <- estimate_no_covariates(data, n_lists = 2)

    expect_named(result, c("listpair", "n", "sigma_n", "ci_l", "ci_u"))
    expect_equal(result$listpair, "1,2")
    expect_equal(result$n, 4)
    expect_equal(result$sigma_n, 1.68, tolerance = 0.005)
    expect_equal(result$ci_l, 4)
    expect_equal(result$ci_u, 8)
})

test_that("estimate_capture_prob routes data without covariates", {
    data <- data.frame(
        capture_1 = c(1, 0, 1, 1),
        capture_2 = c(0, 1, 1, 1)
    )

    expect_warning(
        result <- estimate_capture_prob(data, n_lists = 2),
        "nfolds is reduced to 1"
    )
    expect_equal(result, estimate_no_covariates(data, n_lists = 2))
})

test_that("estimate_capture_prob routes data with covariates", {
    data <- data.frame(
        capture_1 = c(1, 0, 1, 1),
        capture_2 = c(0, 1, 1, 1),
        covariate = c(0, 1, 0, 1)
    )

    result <- estimate_capture_prob(
        data,
        n_lists = 2,
        nfolds = 1,
        method = "plugin",
        func = "logit",
        margin = 0.01,
        seed = 7
    )

    expect_named(result, c("listpair", "n", "sigma_n", "ci_l", "ci_u"))
    expect_equal(result$listpair, "1,2")
    expect_equal(result$n, 4)
})

test_that("estimate_capture_prob validates binary capture columns", {
    data <- data.frame(
        capture_1 = c(1, 0, 2),
        capture_2 = c(0, 1, 1)
    )

    expect_error(
        estimate_capture_prob(data, n_lists = 2),
        "validate_binary_cols"
    )
})

test_that("estimate_capture_prob rejects empty data", {
    data <- data.frame(capture_1 = integer(), capture_2 = integer())

    expect_warning(
        expect_error(
            estimate_capture_prob(data, n_lists = 2),
            "Data must have at least one row"
        ),
        "nfolds is reduced to 1"
    )
})

test_that("estimate_psi returns plugin estimate and variance", {
    nuisance_functions <- list(
        q_1 = c(0.6, 0.7),
        q_2 = c(0.7, 0.8),
        q_12 = c(0.2, 0.3)
    )

    result <- crcsim:::estimate_psi(
        method = "plugin",
        nuisance_functions = nuisance_functions,
        y_j = c(0, 1),
        y_k = c(1, 1),
        margin = 0.1,
        n_lists = 3,
        n = 2
    )

    expect_equal(result$capture_probability, 1.6)
    expect_equal(result$variance, 2.777778, tolerance = 0.000001)
})

test_that("estimate_psi estimates psi with TMLE", {
    nuisance_functions <- list(
        q_1 = c(0.55, 0.60, 0.50, 0.65, 0.45, 0.70, 0.58, 0.62),
        q_2 = c(0.60, 0.50, 0.65, 0.55, 0.70, 0.48, 0.62, 0.57),
        q_12 = c(0.25, 0.20, 0.30, 0.22, 0.28, 0.18, 0.24, 0.26)
    )

    result <- suppressWarnings(crcsim:::estimate_psi(
        method = "tmle",
        nuisance_functions = nuisance_functions,
        y_j = c(0, 1, 0, 1, 1, 0, 1, 0),
        y_k = c(1, 0, 1, 1, 0, 1, 0, 1),
        margin = 0.05,
        n_lists = 3,
        n = 8
    ))

    expect_named(result, c("capture_probability", "variance"))
    expect_true(is.finite(result$capture_probability))
    expect_true(is.finite(result$variance))
    expect_equal(result$capture_probability, 2.357701, tolerance = 0.000001)
    expect_equal(result$variance, 3.41156, tolerance = 0.00001)
})

test_that("estimate_psi uses TMLE and doubly robust methods", {
    nuisance_functions <- list(
        q_1 = c(0.6, 0.7),
        q_2 = c(0.7, 0.8),
        q_12 = c(0.2, 0.3)
    )

    result <- crcsim:::estimate_psi(
        method = "doubly_robust",
        nuisance_functions = nuisance_functions,
        y_j = c(0, 1),
        y_k = c(1, 1),
        margin = 0.1,
        n_lists = 3,
        n = 2
    )

    expect_equal(
        result$capture_probability,
        1.333333,
        tolerance = 0.000001
    )
    expect_equal(result$variance, 2.777778, tolerance = 0.000001)
})