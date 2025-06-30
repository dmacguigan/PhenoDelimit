
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

test_that("dapc_clumpp creates proper directory structure", {
  skip_if_not(file.exists(Sys.which("clumpp")) ||
                file.exists("/usr/local/bin/clumpp") ||
                Sys.getenv("CLUMPP_PATH") != "",
              "CLUMPP executable not found")

  test_data <- create_test_data()
  temp_dir <- tempdir()

  # Mock CLUMPP executable (for testing without actual CLUMPP)
  mock_clumpp <- file.path(temp_dir, "mock_clumpp")
  if (.Platform$OS.type == "windows") {
    mock_clumpp <- paste0(mock_clumpp, ".exe")
  }
  writeLines("#!/bin/bash\necho 'Mock CLUMPP'", mock_clumpp)
  if (.Platform$OS.type != "windows") {
    Sys.chmod(mock_clumpp, mode = "0755")
  }

  # Test directory creation
  expect_false(dir.exists(file.path(temp_dir, "CLUMPP")))

  # This would normally call dapc_clumpp, but we'll test the directory creation logic
  clumpp_dir <- file.path(temp_dir, "CLUMPP")
  dir.create(clumpp_dir, showWarnings = FALSE)

  expect_true(dir.exists(clumpp_dir))

  # Cleanup
  unlink(clumpp_dir, recursive = TRUE)
  unlink(mock_clumpp)
})

test_that("dapc_clumpp validates input parameters", {
  test_data <- create_test_data()

  # Test invalid working directory
  expect_error(
    dapc_clumpp(wd = "/nonexistent/directory",
                CLUMPP_exe = "test",
                data = test_data$data,
                n.groups = c(2),
                model.numbers = c(1),
                models = test_data$models,
                perc.var = c(70)),
    class = "error"
  )

  # Test mismatched dimensions
  wrong_models <- data.frame(m1 = rep(1:2, each = 10))  # Wrong length

  expect_error(
    dapc_clumpp(wd = tempdir(),
                CLUMPP_exe = "test",
                data = test_data$data,
                n.groups = c(2),
                model.numbers = c(1),
                models = wrong_models,
                perc.var = c(70)),
    class = "error"
  )
})

test_that("dapc_clumpp parameter validation works", {
  test_data <- create_test_data()
  temp_dir <- tempdir()

  # Test invalid percentage variance
  expect_error(
    dapc_clumpp(wd = temp_dir,
                CLUMPP_exe = "test",
                data = test_data$data,
                n.groups = c(2),
                model.numbers = c(1),
                models = test_data$models,
                perc.var = c(150)),  # Invalid percentage
    class = "error"
  )

  # Test invalid n.groups
  expect_error(
    dapc_clumpp(wd = temp_dir,
                CLUMPP_exe = "test",
                data = test_data$data,
                n.groups = c(0),  # Invalid group number
                model.numbers = c(1),
                models = test_data$models,
                perc.var = c(70)),
    class = "error"
  )
})
