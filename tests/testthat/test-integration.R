test_that("complete workflow can be tested with mock data", {
  skip_if_not(file.exists(Sys.which("CLUMPP")),
              "Integration test - requires CLUMPP executable")

  library(reshape2)
  library(RColorBrewer)

  # This test would run the complete workflow but is skipped
  # unless CLUMPP is available
  temp_dir <- tempdir()

  create_test_data <- function(wd){
    setwd(wd)

    # how many clusters to generate?
    true_clust = 4
    # how many samples total?
    samples = 300
    # how many variables to generate?
    nVar = 20
    # small variable sd values
    var_sd_min_small = 1
    var_sd_max_small = 3
    # large variable sd values
    var_sd_min_large = 3
    var_sd_max_large = 10
    # large range of variable mean values
    var_mean_min_large = 5
    var_mean_max_large = 30
    # small range of variable mean values
    var_mean_min_small = 18
    var_mean_max_small = 22

    # generate variable normal distribution means
    temp_small <- list()
    temp_large <- list()
    for(i in 1:true_clust){
      temp_large[[i]] <- runif(10, var_mean_min_large, var_mean_max_large) # choose random value from uniform distribution between var_mean_min and var_mean_max
    }
    for(i in 1:true_clust){
      temp_small[[i]] <- runif(10, var_mean_min_small, var_mean_max_small) # choose random value from uniform distribution between var_mean_min and var_mean_max
    }
    var_means_small <- as.data.frame(do.call(rbind, temp_small))
    var_means_large <- as.data.frame(do.call(rbind, temp_large))

    var_means <- cbind(var_means_large, var_means_small)

    # generate variable normal distribution standard deviations, use same sd for each variable
    var_sds <- c(runif(10, var_sd_min_small, var_sd_max_small), runif(10, var_sd_min_large, var_sd_max_large))  # choose random value from uniform distribution between var_sd_min and var_sd_max

    # generate draws from normal distribution for each sample
    sim_data <- data.frame()
    # for each sample
    for(i in 1:samples){
      sim_clust = sample(1:true_clust, 1) # choose a random cluster
      sim_data[i,1] <- sim_clust
      for(j in 1:nVar){
        sim_data[i,j+1] <-rnorm(1, mean=var_means[sim_clust,j], sd=var_sds[j]) # generate draw from normal distribution of approporate cluster
      }
    }
    colnames(sim_data) <- c("cluster", paste("var", c(1:nVar), sep=""))

    # write simulated data to file
    write.table(sim_data[,2:(nVar+1)], file = "sim_data.txt", quote=FALSE, row.names = FALSE)

    # create a few "delimitation models" to test
    # merge clusters 1 + 2 (3 groups total)
    m3 <- replace(sim_data$cluster, sim_data$cluster == 1, 2)
    # merge clusters 1 + 2 and 3 + 4 (2 groups total)
    m2 <- replace(m3, m3==4, 3)
    # randomly split cluster 4 into 2 clusters (5 groups total)
    m4 <- sim_data$cluster
    m4[sample(which(sim_data$cluster == 4), length(which(sim_data$cluster == 4))/2)] <- 5
    # randomly take samples from clusters 3 and 4 and create a new cluster (5 groups total)
    m5 <- sim_data$cluster
    matches <- which(sim_data$cluster == 4 | sim_data$cluster == 3)
    m5[sample(matches, length(matches)/2)] <- 5
    # keep 4 clusters but randomize cluster assignment (4 groups total)
    m6 = sample(sim_data$cluster)

    # combine all models and write to file
    models <- as.data.frame(cbind(m1=sim_data$cluster, m2=m2, m3=m3, m4=m4, m5=m5, m6=m6))
    write.table(models, file = "sim_models.txt", quote=FALSE, row.names = FALSE)
  }


  create_test_data(wd = temp_dir)
  expect_true(file.exists(file.path(temp_dir, "sim_data.txt")))
  expect_true(file.exists(file.path(temp_dir, "sim_models.txt")))

  # step 1: K-means clustering, discriminant analysis of principal components, and run CLUMPP
  wd = temp_dir # results will be written to new subdirectory "CLUMPP"
  CLUMPP_exe = "CLUMPP"
  data = read.table("sim_data.txt", header=TRUE)
  n.groups = c(4,2,3,5,5,4)
  model.numbers = c(1:6)
  models = read.table("sim_models.txt", header=TRUE)
  perc.var = c(70,80,90)
  scale = TRUE
  center = TRUE
  dapc_clumpp(wd=wd,
              CLUMPP_exe=CLUMPP_exe,
              data=data,
              n.groups=n.groups,
              model.numbers=model.numbers,
              models=models,
              perc.var=perc.var,
              scale=scale,
              center=center,
              apriori=FALSE)
  expect_true(dir.exists(file.path(temp_dir, "CLUMPP")))
  
  # let's run this step again but using a priori population assignments
  # we'll come back to these results later
  dapc_clumpp(wd=wd,
              CLUMPP_exe=CLUMPP_exe,
              data=data,
              n.groups=n.groups,
              model.numbers=model.numbers,
              models=models,
              perc.var=perc.var,
              scale=scale,
              center=center,
              apriori=TRUE)

  # step 2: summarize CLUMPP
  
  clumpp_results <- read_clumpp_results(wd=wd,
                                      perc.var=perc.var,
                                      model.numbers=model.numbers)
  expect_true(is.data.frame(clumpp_results))

  # step 3: plot H' values to compare delimitation models
  
  colors = c("gray70", "gray30", "black")
  plot.type = "png"
  plot.name = "H_plot_example"
  plot.width = 8
  plot.height = 4

  p <- plot_clumpp_results(wd=wd,
                      clumpp.data=clumpp_results,
                      colors=colors,
                      plot.name = plot.name,
                      plot.type=plot.type,
                      plot.width=plot.width,
                      plot.height=plot.height)
  expect_true(exists("p"))
    
  # step 4: permutation test of significance for H' values
  
  n.groups = c(4,2,3,5,5,4)
  model.numbers = c(1:6)
  scale = TRUE
  center = TRUE
  dapc_clumpp_permuted(wd=wd,
                     CLUMPP_exe=CLUMPP_exe,
                     data=data,
                     n.groups=n.groups,
                     model.numbers=model.numbers,
                     models=models,
                     perc.var=perc.var,
                     permutations=10,
                     scale=scale,
                     center=center)
  expect_true(dir.exists(file.path(temp_dir, "CLUMPP_permuted")))

  
  clumpp_perm_df <- read_clumpp_results_permuted(wd=wd,
                                               perc.var=perc.var,
                                               model.numbers=model.numbers,
                                               permutations=10)
  expect_true(is.data.frame(clumpp_perm_df))

  
  p <- H_permutation_plot(wd=wd,
                   clumpp.data=clumpp_results,
                   clumpp.data.permuted=clumpp_perm_df,
                   model.numbers=model.numbers,
                   best.perc.var=70,
                   plot.type="png",
                   plot.prefix="example",
                   plot.width=6,
                   plot.height=4)
  expect_true(exists("p"))
  
  # step 5: bar plots of discriminant analysis assignment probabilities
  
  clumpp.wd = "./CLUMPP"
  sample.names = (1:nrow(models))
  sample.plot.groups = models$m1
  sample.plot.groups.order = c(1,2,3,4)
  sample.order = (nrow(models):1)
  best.perc.var = 90
  best.model.number = 1
  plot.type = "png"
  plot.width = 20
  plot.height = 4
  colors = brewer.pal(n = 5, name = "Set1")
  border.color = "gray"

  p <- assign_probs_barplot(wd = wd,
                      clumpp.wd = clumpp.wd,
                      sample.names = sample.names,
                      sample.plot.groups = sample.plot.groups,
                      sample.plot.groups.order = sample.plot.groups.order,
                      #sample.order = sample.order,
                      best.perc.var = best.perc.var,
                      best.model.number = best.model.number,
                      plot.type = plot.type,
                      plot.width = plot.width,
                      plot.height = plot.height,
                      colors = colors,
                      border.color = border.color,
                      apriori=FALSE)
  expect_true(exists("p"))

  # we can also creat a bar plot based on results using a priori assignment of individuals to clusters
  p <- assign_probs_barplot(wd = wd,
                      clumpp.wd = clumpp.wd,
                      sample.names = sample.names,
                      sample.plot.groups = sample.plot.groups,
                      sample.plot.groups.order = sample.plot.groups.order,
                      #sample.order = sample.order,
                      best.perc.var = best.perc.var,
                      best.model.number = best.model.number,
                      plot.type = plot.type,
                      plot.width = plot.width,
                      plot.height = plot.height,
                      colors = colors,
                      border.color = border.color,
                      apriori=TRUE)
  expect_true(exists("p"))
  
  # step 6: scatter plot or density plot of discriminant axes

  plot.width = 6
  plot.height = 6
  colors = brewer.pal(n = 5, name = "Set1")
  border.color = "gray"
  shapes=c(1:4)
  x.axis=1
  y.axis=2

  p <- plot_discriminant_axes(wd = wd,
                        clumpp.wd = clumpp.wd,
                        sample.plot.groups = sample.plot.groups,
                        sample.plot.groups.order = sample.plot.groups.order,
                        best.perc.var = best.perc.var,
                        best.model.number = best.model.number,
                        plot.type = plot.type,
                        plot.width = plot.width,
                        plot.height = plot.height,
                        colors = colors,
                        shapes = shapes,
                        x.axis = x.axis,
                        y.axis = y.axis,
                        apriori=FALSE)
  expect_true(exists("p"))

y.axis=3
p <- plot_discriminant_axes(wd = wd,
                       clumpp.wd = clumpp.wd,
                       sample.plot.groups = sample.plot.groups,
                       sample.plot.groups.order = sample.plot.groups.order,
                       best.perc.var = best.perc.var,
                       best.model.number = best.model.number,
                       plot.type = plot.type,
                       plot.width = plot.width,
                       plot.height = plot.height,
                       colors = colors,
                       shapes = shapes,
                       x.axis = x.axis,
                       y.axis = y.axis,
                       apriori=FALSE)
  expect_true(exists("p"))
  
  x.axis=1
  y.axis=2
  p <- plot_discriminant_axes(wd = wd,
                        clumpp.wd = clumpp.wd,
                        sample.plot.groups = sample.plot.groups,
                        sample.plot.groups.order = sample.plot.groups.order,
                        best.perc.var = best.perc.var,
                        best.model.number = best.model.number,
                        plot.type = plot.type,
                        plot.width = plot.width,
                        plot.height = plot.height,
                        colors = colors,
                        shapes = shapes,
                        x.axis = x.axis,
                        y.axis = y.axis,
                        apriori=TRUE)
  expect_true(exists("p"))
  
  y.axis=3
  p <- plot_discriminant_axes(wd = wd,
                        clumpp.wd = clumpp.wd,
                        sample.plot.groups = sample.plot.groups,
                        sample.plot.groups.order = sample.plot.groups.order,
                        best.perc.var = best.perc.var,
                        best.model.number = best.model.number,
                        plot.type = plot.type,
                        plot.width = plot.width,
                        plot.height = plot.height,
                        colors = colors,
                        shapes = shapes,
                        x.axis = x.axis,
                        y.axis = y.axis,
                        apriori=TRUE)
  expect_true(exists("p"))
  
  # step 7: plot discriminant axis loadings and write table of variable contributions and loadings
  plot.width = 6
  plot.height = 6

  # Variable loading on discriminant axis 1
  axis=1

  p <- discriminant_loading(wd = wd, clumpp.wd = clumpp.wd,
                      best.perc.var = best.perc.var, best.model.number = best.model.number,
                      plot.type = plot.type, plot.width = plot.width, plot.height = plot.height, axis = axis)
  expect_true(exists("p"))
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
