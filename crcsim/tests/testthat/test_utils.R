test_that("build_listpairs returns unique upper-triangular list pairs", {
    result <- build_listpairs(4)

    expect_type(result, "character")
    expect_length(result, 6)
    expect_identical(
        result,
        c("1,2", "1,3", "1,4", "2,3", "2,4", "3,4")
    )
})

test_that("build_listpairs returns one pair for two lists", {
    expect_identical(build_listpairs(2), "1,2")
})

test_that("build_listpairs returns no pairs for one list", {
    expect_length(build_listpairs(1), 0)
})

test_that("test_list_overlap returns NULL for sufficient overlap", {
    result <- test_list_overlap(
        c(1, 1, 0, 0),
        c(1, 0, 1, 0),
        margin = 0.2
    )

    expect_null(result)
})

test_that("test_list_overlap warns for insufficient overlap", {
    expect_warning(
        result <- test_list_overlap(
            c(1, 0, 0, 0),
            c(1, 0, 0, 0),
            margin = 0.3
        ),
        "less than 0.3"
    )
    expect_null(result)
})

test_that("validate_binary_cols accepts binary columns", {
    data <- data.frame(
        capture_1 = c(0, 1, 0),
        capture_2 = c(1, 1, 0),
        value = c(2, 3, 4)
    )

    expect_true(crcsim:::validate_binary_cols(data, end = 2))
    expect_true(crcsim:::validate_binary_cols(data, start = 2, end = 2))
})

test_that("validate_binary_cols rejects non-binary columns", {
    data <- data.frame(
        capture_1 = c(0, 1, 2),
        capture_2 = c(1, 1, 0),
        value = c(2, 3, 4)
    )

    expect_false(crcsim:::validate_binary_cols(data, end = 2))
    expect_false(crcsim:::validate_binary_cols(data, start = 3, end = 3))
})

test_that("validate_binary_cols validates only the requested column range", {
    data <- data.frame(
        metadata = c(2, 3),
        capture_1 = c(0, 1),
        capture_2 = c(1, 0)
    )

    expect_true(crcsim:::validate_binary_cols(data, start = 2, end = 3))
})

test_that("zero_matrix creates a zero-filled named matrix", {
    result <- zero_matrix(
        nrow = 2,
        cols = c("first", "second", "third"),
        rownames = c("a", "b")
    )

    expect_true(is.matrix(result))
    expect_equal(dim(result), c(2, 3))
    expect_true(all(result == 0))
    expect_identical(rownames(result), c("a", "b"))
    expect_identical(colnames(result), c("first", "second", "third"))
})

test_that("zero_matrix supports omitted row names", {
    result <- zero_matrix(3, c("x", "y"))

    expect_equal(dim(result), c(3, 2))
    expect_null(rownames(result))
    expect_identical(colnames(result), c("x", "y"))
    expect_true(all(result == 0))
})

test_that("zero_matrix supports zero columns", {
    result <- zero_matrix(2, character(), rownames = c("a", "b"))

    expect_equal(dim(result), c(2, 0))
    expect_identical(rownames(result), c("a", "b"))
    expect_length(colnames(result), 0)
})
