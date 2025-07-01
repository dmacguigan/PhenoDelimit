#' Run DAPC on permuted data
#' 
#' Permutes the data, then performs performs K-means clustering, runs discriminant analysis, and generates "indfiles" and "paramfiles" for CLUMPP based on user supplied info.
#'
#' @param wd working directory to store CLUMPP files, results will be written to new subdirectory "CLUMPP_permuted"
#' @param CLUMPP_exe full path to CLUMPP executable
#' @param data data frame containing phyotypic data, one row per individual, one column per trait, no missing values
#' @param n.groups vector of the number of groups (populations, species, etc) in each delimitation model
#' @param model.numbers vector containing delimitaiton model numbers
#' @param models data frame containing species delimitation models, one row per individual, one column per model
#' @param perc.var vector containing cumulative percentages of variance to retain for discriminant analyses
#' @param permutations number of permutations to perform
#' @param scale scale PCA? Highly recommended unless you transform data prior to analysis. TRUE or FALSE
#' @param center center PCA? Highly recommended unless you transform data prior to analysis. TRUE or FALSE
#' 
#' @export

dapc_clumpp_permuted <- function(wd, CLUMPP_exe, data, n.groups, model.numbers, models, perc.var, permutations, scale=TRUE, center=TRUE){
  # input error handling
	if(any(perc.var <= 0) || (any(perc.var > 100))){
		stop("perc.var values must be greater than 0 and less than 100")
	}
	if(any(n.groups <= 0)){
		stop("n.groups values must be greater than 0")
	}
  if(permutations <=0 ){
		stop("number of permutations must be greater than 0")
	}

  dir.create(file.path(wd, "CLUMPP_permuted"), showWarnings = FALSE)
  setwd(paste0(wd, "/CLUMPP_permuted"))
  for(i in 1:length(n.groups)){
    for(j in 1:length(perc.var)){
      for(p in 1:permutations){
        p_data <- data[sample(nrow(data)),] # permute the data by shuffling rows
        run_dapc_clumpp_permutation(data=p_data, CLUMPP_exe=CLUMPP_exe, modelNumber=model.numbers[i], k=n.groups[i], model=as.factor(models[,i]), perc.var=perc.var[j], center=center, scale=scale, permutation=p)
      }
    }
  }
}

#' Run DAPC CLUMPP permutation with arguments passed from dapc_clumpp_permuted.
#'
#' @param data data
#' @param CLUMPP_exe CLUMPP_exe
#' @param modelNumber modelNumber
#' @param k k
#' @param model model
#' @param perc.var perc.var
#' @param scale scale
#' @param center center
#' @param permutation permutation
#' 
#' @noRd
# 

run_dapc_clumpp_permutation <- function(data, CLUMPP_exe, modelNumber, k, model, perc.var, scale, center, permutation){
  # make sure model variable is OK format
  model <- as.numeric(as.factor(model))
  # K-means DAPC
  kmeans.cluster <- adegenet::find.clusters(data, max.n.clust = 30, n.start = 1000, n.iter=1e6, n.pca=10000, n.clust=k, center=center, scale=scale) #retain all PCs
  dapc.kmeans <- adegenet::dapc(data, kmeans.cluster$grp, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.var, n.da=10000, center=center, scale=scale)

  # write K-means group assignment to file
  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_Kmeans.grp.p", permutation, ".txt", sep=""), x = dapc.kmeans$grp,
              row.names = FALSE, col.names = FALSE)

  # write DAPC variable contributions to file
  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_var.contr.p", permutation, ".txt", sep=""), x = dapc.kmeans$var.contr)

  # write DAPC variable loadings to file
  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_var.load.p", permutation, ".txt", sep=""), x = dapc.kmeans$var.load)

  # write DAPC individual coordinates to file
  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_ind.coord.p", permutation, ".txt", sep=""), x = dapc.kmeans$ind.coord,
              row.names=FALSE)

  # write DAPC assignment probabilities to file
  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_posteriors.p", permutation, ".txt", sep=""), x = dapc.kmeans$posterior,
              row.names = FALSE, col.names = FALSE)

  # pop
  pop.model <- list()
  for(i in 1:k){
    pop.model[[i]] <-  ifelse(model == i, 1, 0)
  }
  pop.assign <- as.data.frame(do.call(cbind, pop.model))

  # create CLUMPP indfile
  d1 <- data.frame(1:nrow(pop.assign),
                   1:nrow(pop.assign),
                   rep("(x)", nrow(pop.assign)),
                   model,
                   rep(":", nrow(pop.assign)),
                   pop.assign)
  d2 <- data.frame(1:nrow(dapc.kmeans$posterior),
                   1:nrow(dapc.kmeans$posterior),
                   rep("(x)", nrow(dapc.kmeans$posterior)),
                   model,
                   rep(":", nrow(dapc.kmeans$posterior)),
                   dapc.kmeans$posterior)
  colnames(d2) <- colnames(d1)
  d <- rbind(d1, d2)
  # create CLUMPP indfile
  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, ".p", permutation, ".indfile", sep=""), x=d,
              quote=FALSE, col.names = FALSE, row.names=FALSE)

  # create CLUMPP paramfile
  f <- paste("m", modelNumber, "_perVar-", perc.var, ".p", permutation, ".paramfile", sep="")
  file <- file(f, "wb")
  cat("DATATYPE 0", file=file, sep="\n")
  cat(paste("INDFILE m", modelNumber, "_perVar-", perc.var, ".p", permutation, ".indfile", sep=""), file=file, append=TRUE, sep="\n")
  cat(paste("OUTFILE m", modelNumber, "_perVar-", perc.var, ".p", permutation, ".outfile", sep=""), file=file, append=TRUE, sep="\n")
  cat(paste("MISCFILE m", modelNumber, "_perVar-", perc.var, ".p", permutation, ".miscfile", sep=""), file=file, append=TRUE, sep="\n")
  cat(paste("K ", ncol(pop.assign), sep=""), file=file, append=TRUE, sep="\n")
  cat(paste("C ", nrow(pop.assign), sep=""), file=file, append=TRUE, sep="\n")
  cat("R 2\nM 1\nW 0\nS 2\nPRINT_PERMUTED_DATA 0\nPRINT_EVERY_PERM 0\nPRINT_RANDOM_INPUTORDER 0\nOVERRIDE_WARNINGS 0\nORDER_BY_RUN 1", file=file, append=TRUE, sep="\n")

  system(paste(CLUMPP_exe, f, sep=" "), ignore.stdout = TRUE, ignore.stderr = TRUE)

  close(file)
}
