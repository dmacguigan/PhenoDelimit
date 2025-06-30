
test_that("complete workflow can be tested with mock data", {
  skip("Integration test - requires CLUMPP executable")

  # This test would run the complete workflow but is skipped
  # unless CLUMPP is available
  temp_dir <- tempdir()
  test_data <- create_test_data()

  # Full workflow test would go here
  # 1. Run dapc_clumpp
  # 2. Read results
  # 3. Create plots
  # 4. Run permutations
  # 5. Create permutation plots
})

test_that("file I/O operations work correctly", {
  temp_dir <- tempdir()

  # Test writing and reading CSV files
  test_df <- data.frame(
    model = c(1, 2, 3),
    H_prime = c(0.8, 0.7, 0.9),
    perc_var = c(70, 80, 90)
  )

  test_file <- file.path(temp_dir, "test_output.csv")
  write.csv(test_df, test_file, row.names = FALSE)

  expect_true(file.exists(test_file))

  read_df <- read.csv(test_file)
  expect_equal(nrow(read_df), 3)
  expect_equal(names(read_df), names(test_df))

  # Cleanup
  unlink(test_file)
})

test_that("directory operations work correctly", {
  temp_dir <- tempdir()
  test_subdir <- file.path(temp_dir, "test_clumpp")

  # Test directory creation
  expect_false(dir.exists(test_subdir))
  dir.create(test_subdir, showWarnings = FALSE)
  expect_true(dir.exists(test_subdir))

  # Test nested directory creation
  nested_dir <- file.path(test_subdir, "nested")
  dir.create(nested_dir, showWarnings = FALSE)
  expect_true(dir.exists(nested_dir))

  # Cleanup
  unlink(test_subdir, recursive = TRUE)
  expect_false(dir.exists(test_subdir))
})
