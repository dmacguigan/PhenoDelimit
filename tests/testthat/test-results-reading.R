# tests/testthat/test-results-reading.R

test_that("read_clumpp_results handles missing files gracefully", {
  temp_dir <- tempdir()

  # Test with non-existent directory
  expect_error(
    read_clumpp_results(wd = file.path(temp_dir, "nonexistent"),
                        perc.var = c(70, 80),
                        model.numbers = c(1, 2)),
    class = "error"
  )
})

test_that("read_clumpp_results validates parameters", {
  temp_dir <- tempdir()

  # Test invalid percentage variance
  expect_error(
    read_clumpp_results(wd = temp_dir,
                        perc.var = c(-10, 150),  # Invalid values
                        model.numbers = c(1, 2)),
    class = "error"
  )

  # Test empty model numbers
  expect_error(
    read_clumpp_results(wd = temp_dir,
                        perc.var = c(70),
                        model.numbers = c()),  # Empty vector
    class = "error"
  )
})

test_that("read_clumpp_results_permuted validates parameters", {
  temp_dir <- tempdir()

  # Test invalid permutation number
  expect_error(
    read_clumpp_results_permuted(wd = temp_dir,
                                 perc.var = c(70),
                                 model.numbers = c(1),
                                 permutations = 0),  # Invalid
    class = "error"
  )

  expect_error(
    read_clumpp_results_permuted(wd = temp_dir,
                                 perc.var = c(70),
                                 model.numbers = c(1),
                                 permutations = -5),  # Invalid
    class = "error"
  )
})

# Mock CLUMPP results for testing
create_mock_clumpp_results <- function(temp_dir) {
  clumpp_dir <- file.path(temp_dir, "CLUMPP")
  dir.create(clumpp_dir, showWarnings = FALSE)

  # Create mock result files
  mock_result <- data.frame(
    H_prime = c(0.85, 0.92, 0.78),
    model = c(1, 1, 1),
    perc_var = c(70, 80, 90)
  )

  write.csv(mock_result,
            file.path(clumpp_dir, "H_results_model1.csv"),
            row.names = FALSE)

  return(clumpp_dir)
}

test_that("mock CLUMPP results can be read", {
  temp_dir <- tempdir()
  clumpp_dir <- create_mock_clumpp_results(temp_dir)

  expect_true(dir.exists(clumpp_dir))
  expect_true(file.exists(file.path(clumpp_dir, "H_results_model1.csv")))

  # Test reading the mock file
  mock_data <- read.csv(file.path(clumpp_dir, "H_results_model1.csv"))
  expect_equal(nrow(mock_data), 3)
  expect_true("H_prime" %in% names(mock_data))

  # Cleanup
  unlink(clumpp_dir, recursive = TRUE)
})
