# Create mock CLUMPP data for plotting tests
create_mock_plot_data <- function() {
  data.frame(
    model = rep(1:3, each = 3),
    perc.var = rep(c(70, 80, 90), 3),
    H_prime = c(0.85, 0.90, 0.95, 0.70, 0.75, 0.80, 0.60, 0.65, 0.70)
  )
}

test_that("plot_clumpp_results validates input data", {
  temp_dir <- tempdir()

  # Test with invalid data
  invalid_data <- data.frame(x = 1:5, y = 6:10)  # Wrong structure

  expect_error(
    plot_clumpp_results(wd = temp_dir,
                        clumpp.data = invalid_data,
                        colors = c("red", "blue"),
                        plot.name = "test",
                        plot.type = "png"),
    class = "error"
  )
})

test_that("plot_clumpp_results validates parameters", {
  temp_dir <- tempdir()
  mock_data <- create_mock_plot_data()

  # Test invalid plot type
  expect_error(
    plot_clumpp_results(wd = temp_dir,
                        clumpp.data = mock_data,
                        colors = c("red", "blue"),
                        plot.name = "test",
                        plot.type = "invalid"),  # Invalid type
    class = "error"
  )

  # Test invalid dimensions
  expect_error(
    plot_clumpp_results(wd = temp_dir,
                        clumpp.data = mock_data,
                        colors = c("red", "blue"),
                        plot.name = "test",
                        plot.type = "png",
                        plot.width = -5),  # Invalid width
    class = "error"
  )
})

test_that("assign_probs_barplot validates parameters", {
  temp_dir <- tempdir()
  clumpp_dir <- file.path(temp_dir, "CLUMPP")
  dir.create(clumpp_dir, showWarnings = FALSE)

  # Test invalid sample names
  expect_error(
    assign_probs_barplot(wd = temp_dir,
                         clumpp.wd = clumpp_dir,
                         sample.names = c(),  # Empty
                         sample.plot.groups = c(1, 2, 3),
                         best.perc.var = 70,
                         best.model.number = 1),
    class = "error"
  )

  # Test mismatched lengths
  expect_error(
    assign_probs_barplot(wd = temp_dir,
                         clumpp.wd = clumpp_dir,
                         sample.names = 1:10,
                         sample.plot.groups = c(1, 2, 3),  # Wrong length
                         best.perc.var = 70,
                         best.model.number = 1),
    class = "error"
  )

  # Cleanup
  unlink(clumpp_dir, recursive = TRUE)
})

test_that("plot_discriminant_axes validates parameters", {
  temp_dir <- tempdir()
  clumpp_dir <- file.path(temp_dir, "CLUMPP")

  # Test invalid axis specification
  expect_error(
    plot_discriminant_axes(wd = temp_dir,
                           clumpp.wd = clumpp_dir,
                           sample.plot.groups = c(1, 2, 3, 4),
                           best.perc.var = 70,
                           best.model.number = 1,
                           x.axis = 0),  # Invalid axis
    class = "error"
  )

  expect_error(
    plot_discriminant_axes(wd = temp_dir,
                           clumpp.wd = clumpp_dir,
                           sample.plot.groups = c(1, 2, 3, 4),
                           best.perc.var = 70,
                           best.model.number = 1,
                           x.axis = 1,
                           y.axis = 1),  # Same axis for x and y
    class = "error"
  )
})

test_that("discriminant_loading validates parameters", {
  temp_dir <- tempdir()
  clumpp_dir <- file.path(temp_dir, "CLUMPP")

  # Test invalid axis
  expect_error(
    discriminant_loading(wd = temp_dir,
                         clumpp.wd = clumpp_dir,
                         best.perc.var = 70,
                         best.model.number = 1,
                         axis = 0),  # Invalid axis
    class = "error"
  )
})
