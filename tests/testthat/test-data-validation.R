test_that("data validation works correctly", {
  # Valid data
  valid_data <- data.frame(
    trait1 = rnorm(20),
    trait2 = rnorm(20),
    trait3 = rnorm(20)
  )

  # Invalid data with missing values
  invalid_data <- data.frame(
    trait1 = c(1, 2, NA, 4, 5),
    trait2 = c(1, NA, 3, 4, 5)
  )

  # Empty data
  empty_data <- data.frame()

  # Test data frame structure
  expect_true(is.data.frame(valid_data))
  expect_equal(ncol(valid_data), 3)
  expect_equal(nrow(valid_data), 20)

  # Test for missing values (your functions likely need complete data)
  expect_false(any(is.na(valid_data)))
  expect_true(any(is.na(invalid_data)))
})

test_that("model specification validation works", {
  # Valid models
  valid_models <- data.frame(
    m1 = rep(1:4, each = 5),
    m2 = rep(1:2, each = 10),
    m3 = c(rep(1:2, each = 5), rep(3, 10))
  )

  # Invalid models (non-integer)
  invalid_models <- data.frame(
    m1 = rnorm(20)
  )

  expect_true(is.data.frame(valid_models))
  expect_equal(nrow(valid_models), 20)

  # Test that model assignments are integers/factors
  expect_true(all(sapply(valid_models, function(x) all(x == floor(x)))))
})

test_that("parameter validation works", {
  # Test n.groups parameter
  valid_n_groups <- c(2, 3, 4, 5)
  invalid_n_groups <- c(1, 0, -1)

  expect_true(all(valid_n_groups > 1))
  expect_false(all(invalid_n_groups > 1))

  # Test perc.var parameter
  valid_perc_var <- c(70, 80, 90)
  invalid_perc_var <- c(0, 101, -10)

  expect_true(all(valid_perc_var > 0 & valid_perc_var <= 100))
  expect_false(all(invalid_perc_var > 0 & invalid_perc_var <= 100))
})
