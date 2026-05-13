library(crcsim)

set.seed(2026)

data_table <- create_data(
    n_individuals = 10000,
    p_captures = c(0.45, 0.3, 0.2),
    p_strata = 1
)

capture_columns <- c("capture_1", "capture_2", "capture_3")

# crc expects only observed rows, so remove all-zero capture rows
observed_table <- data_table[rowSums(data_table[, capture_columns]) != 0, ]

tmle_result <- crc(
    data = observed_table,
    freq_column = "N_ID",
    binary_variables = capture_columns,
    formula_selection = "tmle",
    opts_tmle = list(
        list_pair = c("capture_1", "capture_2"),
        funcname = "logit",
        nfolds = 2,
        margin = 0.005,
        estimator = "TMLE",
        expansion_mode = "exact",
        sample_size = 5000,
        warn_expanded_rows = 1000000
    )
)

tmle_result
