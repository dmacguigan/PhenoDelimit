
# Helper functions for creating test data across multiple test files

#' Create standardized test dataset
#'
#' @param n_individuals Number of individuals
#' @param n_traits Number of traits
#' @param n_groups Number of groups
#' @return List with data and models
create_standard_test_data <- function(n_individuals = 50,
                                      n_traits = 5,
                                      n_groups = 2) {
  set.seed(42)  # Reproducible

  # Create trait data with group differences
  data <- data.frame(
    lapply(1:n_traits, function(i) {
      group_means <- seq(10, 20, length.out = n_groups)
      trait_values <- numeric(n_individuals)
      individuals_per_group <- n_individuals / n_groups

      for (g in 1:n_groups) {
        start_idx <- (g - 1) * individuals_per_group + 1
        end_idx <- g * individuals_per_group
        trait_values[start_idx:end_idx] <- rnorm(individuals_per_group,
                                                 group_means[g],
                                                 2)
      }
      trait_values
    })
  )

  names(data) <- paste0("trait", 1:n_traits)

  # Create model assignments
  models <- data.frame(
    m1 = rep(1:n_groups, each = n_individuals / n_groups),
    m2 = rep(1, n_individuals)  # Single group model
  )

  list(data = data, models = models)
}

#' Check if CLUMPP is available
#'
#' @return Logical indicating if CLUMPP can be found
clumpp_available <- function() {
  file.exists(Sys.which("clumpp")) ||
    file.exists("/usr/local/bin/clumpp") ||
    Sys.getenv("CLUMPP_PATH") != ""
}

#' Skip test if CLUMPP not available
skip_if_no_clumpp <- function() {
  skip_if_not(clumpp_available(), "CLUMPP executable not found")
}
