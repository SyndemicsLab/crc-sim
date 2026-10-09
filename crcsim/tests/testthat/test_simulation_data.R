test_that("simulate_data returns complete and observed simulation views", {
    result <- simulate_data(
        100,
        3,
        base = homogeneous_capture_spec(rep(0.5, 3)),
        seed = 1
    )
    capture_columns <- result$metadata$capture_columns

    expect_s3_class(result, "crcsim_simulation")
    expect_named(
        result,
        c(
            "full_data",
            "observed_data",
            "unobserved_data",
            "observed_frequency_table",
            "truth",
            "metadata"
        )
    )
    expect_equal(result$truth$n_total, 100L)
    expect_equal(
        result$truth$n_observed + result$truth$n_unobserved,
        result$truth$n_total
    )
    expect_true(all(rowSums(result$observed_data[, capture_columns]) > 0))
    expect_true(all(rowSums(result$unobserved_data[, capture_columns]) == 0))
    expect_equal(
        sum(result$observed_frequency_table$N_ID),
        result$truth$n_observed
    )
})

test_that("simulate_data is reproducible when a seed is supplied", {
    first <- simulate_data(25, 2, seed = 42)
    second <- simulate_data(25, 2, seed = 42)

    expect_identical(first$full_data, second$full_data)
    expect_error(simulate_data(10, 2, seed = 1.5), "integer scalar")
})

test_that("simulate_data supports list-specific covariate capture models", {
    covariate_ranges <- data.frame(
        distribution = c("uniform", "uniform"),
        p1 = c(0, 0),
        p2 = c(1, 1),
        dtype = c("integer", "integer")
    )
    alpha <- c(-1.2, -0.8, -1.5)
    beta <- rbind(
        capture_1 = c(1.0, -0.5),
        capture_2 = c(0.4, 0.8),
        capture_3 = c(-0.7, 1.2)
    )

    result <- simulate_data(
        n_individuals = 100,
        n_captures = 3,
        base = covariate_capture_spec(covariate_ranges, alpha, beta),
        seed = 1
    )

    expect_equal(result$metadata$n_captures, 3L)
    expect_equal(result$truth$alpha, alpha)
    expect_identical(dim(result$truth$beta), c(3L, 2L))
    expect_identical(dim(result$truth$capture_probabilities), c(100L, 3L))
    expect_true(all(result$truth$capture_probabilities > 0))
    expect_true(all(result$truth$capture_probabilities < 1))
    expect_equal(
        result$truth$n_observed + result$truth$n_unobserved,
        result$truth$n_total
    )
})

test_that("simulate_data validates covariate capture model arguments", {
    beta <- matrix(c(0.2, 0.4), nrow = 1)
    ranges <- data.frame(
        distribution = "uniform",
        p1 = 0,
        p2 = 1,
        dtype = "integer"
    )

    expect_error(
        covariate_capture_spec(ranges, 0, beta),
        "one column for each covariate"
    )
    expect_error(
        covariate_capture_spec(
            ranges,
            c(0, 0),
            matrix(c(0.2, 0.4), nrow = 1)
        ),
        "one row"
    )
    expect_error(
        covariate_capture_spec(
            ranges,
            c(0, 0),
            matrix(c(0.2, 0.4), nrow = 1)
        ),
        "one row"
    )
    expect_error(
        covariate_capture_spec(
            ranges,
            0,
            matrix(c(0.2, 0.4), nrow = 1)
        ),
        "one column for each covariate"
    )
})

test_that("simulate_data supports latent heterogeneity", {
    result <- simulate_data(
        n_individuals = 2000,
        n_captures = 3,
        base = homogeneous_capture_spec(rep(0.5, 3)),
        heterogeneity = continuous_heterogeneity_spec(1),
        seed = 6
    )

    expect_equal(result$metadata$scenario, "continuous latent heterogeneity")
    expect_length(result$truth$latent_effects, 2000)
    expect_true(sd(result$truth$latent_effects) > 0.8)
    expect_true(
        result$truth$n_observed + result$truth$n_unobserved ==
            result$truth$n_total
    )
    expect_true(all(result$truth$capture_probabilities > 0))
    expect_true(all(result$truth$capture_probabilities < 1))
})

test_that("simulate_data supports direct pairwise dependence", {
    result <- simulate_data(
        n_individuals = 2000,
        n_captures = 3,
        base = homogeneous_capture_spec(rep(0.5, 3)),
        dependence = pairwise_dependence_spec(c(1, 2), 2),
        seed = 7
    )
    captures <- result$full_data[, c("capture_1", "capture_2")]
    observed_ratio <- mean(captures$capture_1 * captures$capture_2) /
        (mean(captures$capture_1) * mean(captures$capture_2))

    expect_equal(result$metadata$scenario, "direct pairwise dependence")
    expect_equal(result$truth$dependence_pair, c(1, 2))
    expect_gt(observed_ratio, 1.2)
})

test_that("simulate_data supports latent class heterogeneity", {
    class_profiles <- rbind(
        c(0.8, 0.8, 0.8),
        c(0.2, 0.2, 0.2)
    )
    result <- simulate_data(
        n_individuals = 2000,
        n_captures = 3,
        heterogeneity = latent_class_spec(c(3, 1), class_profiles),
        seed = 8
    )
    summary <- summarize_simulation(
        result,
        capture_columns = result$metadata$capture_columns
    )

    expect_equal(result$metadata$scenario, "latent class heterogeneity")
    expect_equal(result$metadata$n_captures, 3L)
    expect_length(result$truth$class_assignments, 2000)
    expect_equal(result$truth$latent_class_probabilities, c(0.75, 0.25))
    expect_equal(result$truth$latent_class_p_captures, class_profiles)
    expect_gt(
        summary$overlaps$observed_to_expected[
            summary$overlaps$scope == "full" &
                summary$overlaps$capture_1 == "capture_1" &
                summary$overlaps$capture_2 == "capture_2"
        ],
        1.08
    )
})

test_that("simulate_data supports graph-structured dependence", {
    edges <- data.frame(
        from = c(1, 2),
        to = c(2, 3),
        strength = c(2, -2)
    )
    result <- simulate_data(
        n_individuals = 2000,
        n_captures = 3,
        base = homogeneous_capture_spec(rep(0.5, 3)),
        dependence = graph_dependence_spec(edges),
        seed = 9
    )
    summary <- summarize_simulation(
        result,
        capture_columns = result$metadata$capture_columns
    )
    overlap <- summary$overlaps[summary$overlaps$scope == "full", ]

    expect_equal(result$metadata$scenario, "graph-structured dependence")
    expect_equal(result$truth$dependence_edges, edges)
    expect_gt(
        overlap$observed_to_expected[
            overlap$capture_1 == "capture_1" &
                overlap$capture_2 == "capture_2"
        ],
        1.2
    )
    expect_lt(
        overlap$observed_to_expected[
            overlap$capture_1 == "capture_2" &
                overlap$capture_2 == "capture_3"
        ],
        0.8
    )
})

test_that("summarize_simulation summarizes all simulation scopes", {
    simulation <- simulate_data(100, 3, seed = 1)
    summary <- summarize_simulation(
        simulation,
        capture_columns = simulation$metadata$capture_columns
    )

    expect_s3_class(summary, "crcsim_summary")
    expect_true(all(
        c(
            "counts",
            "capture_patterns",
            "list_counts",
            "overlaps",
            "covariates",
            "capture_probabilities",
            "metadata"
        ) %in%
            names(summary)
    ))
    expect_setequal(
        summary$counts$scope,
        c("full", "observed", "unobserved")
    )
    expect_equal(
        summary$counts$n_individuals[
            summary$counts$scope == "full"
        ],
        100
    )
    expect_equal(
        sum(summary$capture_patterns$n[
            summary$capture_patterns$scope == "full"
        ]),
        100
    )
    expect_true(all(
        summary$overlaps$jaccard[
            summary$overlaps$scope == "full"
        ] >=
            0
    ))
    expect_equal(summary$metadata$input_type, "crcsim_simulation")
    expect_equal(
        summary$capture_probabilities$mean[
            summary$capture_probabilities$scope == "full"
        ],
        rep(0.5, 3)
    )
})

test_that("summarize_simulation summarizes covariate-dependent probabilities", {
    ranges <- data.frame(
        distribution = c("uniform", "uniform"),
        p1 = c(0, 0),
        p2 = c(1, 1),
        dtype = c("integer", "integer")
    )
    simulation <- simulate_data(
        n_individuals = 40,
        n_captures = 2,
        base = covariate_capture_spec(
            ranges,
            c(-1, 0),
            matrix(c(0.5, -0.25, 0.2, 0.75), nrow = 2)
        ),
        seed = 5
    )
    summary <- summarize_simulation(
        simulation,
        capture_columns = simulation$metadata$capture_columns
    )

    expected <- colMeans(simulation$truth$capture_probabilities)
    expect_equal(
        summary$capture_probabilities$mean[
            summary$capture_probabilities$scope == "full"
        ],
        unname(expected)
    )
    expect_true(all(is.na(summary$capture_probabilities$mean[
        summary$capture_probabilities$scope != "full"
    ])))
})

test_that("summarize_simulation preserves unavailable raw-data scopes", {
    simulation <- simulate_data(100, 3, seed = 2)
    summary <- summarize_simulation(
        simulation$observed_data,
        capture_columns = simulation$metadata$capture_columns,
        scope = "observed"
    )

    expect_equal(summary$metadata$available_scopes, "observed")
    expect_true(all(
        c(
            "full",
            "observed",
            "unobserved"
        ) %in%
            summary$counts$scope
    ))
    expect_true(all(is.na(summary$counts$n_individuals[
        summary$counts$scope != "observed"
    ])))
    expect_true(all(
        summary$capture_patterns$history[
            summary$capture_patterns$scope == "observed"
        ] !=
            "000"
    ))
})

test_that("summarize_simulation includes numeric covariate summaries", {
    ranges <- data.frame(
        distribution = c("uniform", "normal"),
        p1 = c(0, 1),
        p2 = c(1, 0.5),
        dtype = c("integer", "numeric")
    )
    simulation <- simulate_data(
        n_individuals = 50,
        n_captures = 2,
        base = covariate_capture_spec(
            ranges,
            c(0, 0),
            matrix(0, nrow = 2, ncol = 2)
        ),
        seed = 3
    )
    summary <- summarize_simulation(
        simulation,
        capture_columns = c("capture_1", "capture_2")
    )

    expect_setequal(
        summary$covariates$covariate,
        c("covariate_1", "covariate_2")
    )
    expect_true(all(
        summary$covariates$n[
            summary$covariates$scope == "full"
        ] ==
            50
    ))
})

test_that("summarize_simulation handles one capture list", {
    simulation <- simulate_data(20, 1, seed = 4)
    summary <- summarize_simulation(
        simulation,
        capture_columns = "capture_1"
    )

    expect_equal(nrow(summary$overlaps), 3)
    expect_true(all(is.na(summary$overlaps$intersection)))
})

test_that("summarize_simulation summarizes categorical levels and missingness", {
    observed_data <- tibble::tibble(
        capture_1 = c(1, 0, 1, 1),
        capture_2 = c(0, 1, 1, 0),
        sex = factor(c("A", NA, "B", "A"))
    )
    summary <- summarize_simulation(
        observed_data,
        capture_columns = c("capture_1", "capture_2"),
        scope = "observed"
    )

    level_summary <- summary$covariate_levels
    expect_setequal(
        level_summary$level[level_summary$scope == "observed"],
        c("A", "B", "(Missing)")
    )
    expect_equal(
        level_summary$n[
            level_summary$level == "A" &
                level_summary$scope == "observed"
        ],
        2L
    )
    expect_equal(
        level_summary$proportion[
            level_summary$level == "(Missing)" &
                level_summary$scope == "observed"
        ],
        0.25
    )
})
