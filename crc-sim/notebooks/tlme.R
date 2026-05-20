library(crcsim)

set.seed(2026)

data_table <- create_data(
    n_individuals = 10000,
    n_captures = 3,
    n_categories = 2,
    p_captures = c(0.45, 0.3, 0.2),
    p_categories = c(0.4, 0.6)
)

capture_columns <- c("capture_1", "capture_2", "capture_3")

# crc expects only observed rows, so remove all-zero capture rows
observed_table <- data_table[rowSums(data_table[, capture_columns]) != 0, ] |>
    all_int_cols_to_numeric()

qhat <- popsize(
    data = observed_table,
    funcname = c("logit"),
    nfolds = 2,
    margin = 0.005
)

psin_estimate <- popsize(
    data = observed_table,
    getnuis = qhat$nuis,
    idfold = qhat$idfold
)

tmle_result
