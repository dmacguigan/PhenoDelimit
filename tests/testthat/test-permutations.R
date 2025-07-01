
test_that("dapc_clumpp_permuted validates parameters", {
  skip_if_not(file.exists(Sys.which("CLUMPP")) ||
                Sys.getenv("CLUMPP_PATH") != "",
              "CLUMPP executable not found")

  # Create test data
  create_test_data <- function() {
    set.seed(123)
    data <- data.frame(
      trait1 = c(rnorm(25, 10, 2), rnorm(25, 15, 2)),
      trait2 = c(rnorm(25, 5, 1), rnorm(25, 8, 1)),
      trait3 = c(rnorm(25, 20, 3), rnorm(25, 25, 3))
    )

    models <- data.frame(
      m1 = rep(1:2, each = 25),
      m2 = rep(1, 50)  # Single group model
    )

    list(data = data, models = models)
  }

  test_data <- create_test_data()
  temp_dir <- tempdir()

  # Test invalid permutation number
  expect_error(
    dapc_clumpp_permuted(wd = temp_dir,
                         CLUMPP_exe = "CLUMPP",
                         data = test_data$data,
                         n.groups = c(2),
                         model.numbers = c(1),
                         models = test_data$models,
                         perc.var = c(70),
                         permutations = 0),  # Invalid
    class = "error"
  )

  expect_error(
    dapc_clumpp_permuted(wd = temp_dir,
                         CLUMPP_exe = "CLUMPP",
                         data = test_data$data,
                         n.groups = c(2),
                         model.numbers = c(1),
                         models = test_data$models,
                         perc.var = c(70),
                         permutations = -5),  # Invalid
    class = "error"
  )
})

test_that("H_permutation_plot validates input data", {
  temp_dir <- tempdir()

  # Mock observed data
  observed_data <- data.frame(
    model = c(1, 2),
    perc.var = c(70, 70),
    H_prime = c(0.85, 0.75)
  )

  # Mock permuted data
  permuted_data <- data.frame(
    model = rep(c(1, 2), each = 100),
    perc.var = rep(70, 200),
    H_prime = c(rnorm(100, 0.5, 0.1), rnorm(100, 0.4, 0.1)),
    permutation = rep(1:100, 2)
  )

  # Test mismatched data
  wrong_observed <- data.frame(model = 3, perc.var = 70, H_prime = 0.9)

  expect_error(
    H_permutation_plot(wd = temp_dir,
                       clumpp.data = wrong_observed,
                       clumpp.data.permuted = permuted_data,
                       model.numbers = c(1, 2),
                       best.perc.var = 70),
    class = "error"
  )
})

test_that("permutation significance testing logic works", {
  # Create test permutation data
  observed_h <- 0.85
  permuted_h <- rnorm(100, 0.5, 0.1)

  # Calculate p-value (proportion of permuted values >= observed)
  p_value <- sum(permuted_h >= observed_h) / length(permuted_h)

  expect_true(p_value >= 0 && p_value <= 1)
  expect_true(p_value < 0.05)  # Should be significant

  # Test with non-significant case
  observed_h_low <- 0.45
  p_value_low <- sum(permuted_h >= observed_h_low) / length(permuted_h)
  expect_true(p_value_low > 0.05)  # Should not be significant
})
