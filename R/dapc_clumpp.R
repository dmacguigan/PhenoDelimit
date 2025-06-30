#' Run DAPC on data
#' 
#' Permutes the data, then performs performs K-means clustering, runs discriminant analysis, and generates "indfiles" and "paramfiles" for CLUMPP based on user supplied info.
#' Note that PCA with scale=TRUE is a correlation-based PCA, while PCA with scale=FALSE is a covariance-based PCA.
#' See https://aedin.github.io/PCAworkshop/articles/b_PCA.html.
#' These two types of PCA are appropriate for different datasets.
#' For example, covariance-based PCA is appropriate for morphometric data.
#' See https://www.palass.org/publications/newsletter/palaeomath-101/palaeomath-part-21-principal-warps-relative-warps-and-procrustes-pca.
#' While correlation-based PCA is appropriate for datasets with mixed variable types that are on different scales (e.g. meristic data).
#' See https://stats.stackexchange.com/questions/53/pca-on-correlation-or-covariance for more info.
#' "You tend to use the covariance matrix when the variable scales are similar and the correlation matrix when variables are on different scales.
#'
#' @param wd working directory to store CLUMPP files, results will be written to new subdirectory "CLUMPP"
#' @param CLUMPP_exe full path to CLUMPP executable
#' @param data data frame containing phyotypic data, one row per individual, one column per trait, no missing values
#' @param n.groups vector of the number of groups (populations, species, etc) in each delimitation model
#' @param model.numbers vector containing delimitaiton model numbers
#' @param models data frame containing species delimitation models, one row per individual, one column per model
#' @param perc.var vector containing cumulative percentages of variance to retain for discriminant analyses
#' @param scale scale PCA? Highly recommended unless you transform data prior to analysis. TRUE or FALSE
#' @param center center PCA? Highly recommended unless you transform data prior to analysis. TRUE or FALSE
#' @param apriori do you wish to use apriori individual assignment to species/populations/clusters, or assign individuals using k-means clustering? TRUE or FALSE
#'
#' @export
dapc_clumpp <- function(wd, CLUMPP_exe, data, n.groups, model.numbers, models, perc.var, scale=TRUE, center=TRUE, apriori=FALSE){
  clust.method="kmeans"
  dir.create(file.path(wd, "CLUMPP"), showWarnings = FALSE)
  setwd(paste0(wd, "/CLUMPP"))
  for(i in 1:length(n.groups)){
    for(j in 1:length(perc.var)){
      run_dapc_clumpp(data=data, CLUMPP_exe=CLUMPP_exe, modelNumber=model.numbers[i], k=n.groups[i], model=as.factor(models[,i]), perc.var=perc.var[j], center=center, scale=scale, apriori=apriori, clust.method=clust.method)
    }
  }
}

#' Run DAPC CLUMPP with arguments passed from dapc_clumpp
#'
#' @param data data
#' @param CLUMPP_exe CLUMPP_exe
#' @param modelNumber modelNumber
#' @param k k
#' @param model model
#' @param perc.var perc.var
#' @param scale scale
#' @param center center
#' @param clust.method clust.method
#' @param apriori apriori
#'
#' @noRd
#
run_dapc_clumpp <- function(data, CLUMPP_exe, modelNumber, k, model, perc.var, scale, center, clust.method, apriori){
  if(clust.method == "kmeans"){
	  # make sure model variable is OK format
	  model <- as.numeric(as.factor(model))
	  if(apriori == FALSE){
		  # K-means DAPC
		  kmeans.cluster <- find.clusters(data, max.n.clust = 30, n.start = 1000, n.iter=1e6, n.pca=10000, n.clust=k, center=center, scale=scale) #retain all PCs
		  dapc.kmeans <- dapc(data, kmeans.cluster$grp, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.var, n.da=10000, center=center, scale=scale)

		  # write K-means group assignment to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_Kmeans.grp.txt", sep=""), x = dapc.kmeans$grp,
					  row.names = FALSE, col.names = FALSE)

		  # write DAPC variable contributions to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_var.contr.txt", sep=""), x = dapc.kmeans$var.contr)

		  # write DAPC variable loadings to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_var.load.txt", sep=""), x = dapc.kmeans$var.load)

		  # write DAPC individual coordinates to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_ind.coord.txt", sep=""), x = dapc.kmeans$ind.coord,
					  row.names=FALSE)

		  # write DAPC assignment probabilities to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_posteriors.txt", sep=""), x = dapc.kmeans$posterior,
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
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, ".indfile", sep=""), x=d,
					  quote=FALSE, col.names = FALSE, row.names=FALSE)

		  # create CLUMPP paramfile
		  f <- paste("m", modelNumber, "_perVar-", perc.var, ".paramfile", sep="")
		  file <- file(f, "wb")
		  cat("DATATYPE 0", file=file, sep="\n")
		  cat(paste("INDFILE m", modelNumber, "_perVar-", perc.var, ".indfile", sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("OUTFILE m", modelNumber, "_perVar-", perc.var, ".outfile", sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("MISCFILE m", modelNumber, "_perVar-", perc.var, ".miscfile", sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("K ", ncol(pop.assign), sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("C ", nrow(pop.assign), sep=""), file=file, append=TRUE, sep="\n")
		  cat("R 2\nM 1\nW 0\nS 2\nPRINT_PERMUTED_DATA 0\nPRINT_EVERY_PERM 0\nPRINT_RANDOM_INPUTORDER 0\nOVERRIDE_WARNINGS 0\nORDER_BY_RUN 1", file=file, append=TRUE, sep="\n")

		  system(paste(CLUMPP_exe, f, sep=" "))

		  close(file)
	  }

	  else{
		  # DAPC using a priori assignment of individuals
		  dapc.apriori <- dapc(data, model, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.var, n.da=10000, center=center, scale=scale)

		  # write DAPC variable contributions to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_var.contr.apriori.txt", sep=""), x = dapc.apriori$var.contr)

		  # write DAPC variable loadings to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_var.load.apriori.txt", sep=""), x = dapc.apriori$var.load)

		  # write DAPC individual coordinates to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_ind.coord.apriori.txt", sep=""), x = dapc.apriori$ind.coord,
					  row.names=FALSE)

		  # write DAPC assignment probabilities to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_posteriors.apriori.txt", sep=""), x = dapc.apriori$posterior,
					  row.names = FALSE, col.names = FALSE)
	  }
  } else if(clust.method == "randomforest"){
	  # make sure model variable is OK format
	  model <- as.numeric(as.factor(model))
	  if(apriori == FALSE){
		  # first, perform PCA
		  maxRank <- min(dim(data))
		  pcaX <- dudi.pca(data, center = center, scale = scale, scannf = FALSE, nf=maxRank)
		  pcs <- pcaX$li
		  # then perform random forest clustering
		  rf.fit <- randomForest(x = pcs, y = NULL, ntree = 10000, proximity = TRUE, oob.prox = TRUE)
		  hclust.rf <- hclust(as.dist(1-rf.fit$proximity), method = "ward.D2")
		  rf.cluster = cutree(hclust.rf, k=k)
		  dapc.rf <- dapc(data, rf.cluster, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.var, n.da=10000, center=center, scale=scale)

		  # write RF group assignment to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.grp.txt", sep=""), x = rf.cluster,
					  row.names = FALSE, col.names = FALSE)

		  # write DAPC variable contributions to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.var.contr.txt", sep=""), x = dapc.rf$var.contr)

		  # write DAPC variable loadings to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.var.load.txt", sep=""), x = dapc.rf$var.load)

		  # write DAPC individual coordinates to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.ind.coord.txt", sep=""), x = dapc.rf$ind.coord,
					  row.names=FALSE)

		  # write DAPC assignment probabilities to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.posteriors.txt", sep=""), x = dapc.rf$posterior,
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
		  d2 <- data.frame(1:nrow(dapc.rf$posterior),
						   1:nrow(dapc.rf$posterior),
						   rep("(x)", nrow(dapc.rf$posterior)),
						   model,
						   rep(":", nrow(dapc.rf$posterior)),
						   dapc.rf$posterior)
		  colnames(d2) <- colnames(d1)
		  d <- rbind(d1, d2)
		  # create CLUMPP indfile
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, ".RF.indfile", sep=""), x=d,
					  quote=FALSE, col.names = FALSE, row.names=FALSE)

		  # create CLUMPP paramfile
		  f <- paste("m", modelNumber, "_perVar-", perc.var, ".RF.paramfile", sep="")
		  file <- file(f, "wb")
		  cat("DATATYPE 0", file=file, sep="\n")
		  cat(paste("INDFILE m", modelNumber, "_perVar-", perc.var, ".RF.indfile", sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("OUTFILE m", modelNumber, "_perVar-", perc.var, ".RF.outfile", sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("MISCFILE m", modelNumber, "_perVar-", perc.var, ".RF.miscfile", sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("K ", ncol(pop.assign), sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("C ", nrow(pop.assign), sep=""), file=file, append=TRUE, sep="\n")
		  cat("R 2\nM 1\nW 0\nS 2\nPRINT_PERMUTED_DATA 0\nPRINT_EVERY_PERM 0\nPRINT_RANDOM_INPUTORDER 0\nOVERRIDE_WARNINGS 0\nORDER_BY_RUN 1", file=file, append=TRUE, sep="\n")

		  system(paste("CLUMPP", f, sep=" "))

		  close(file)
	  }

	  else{
		  # first, perform PCA
		  maxRank <- min(dim(data))
		  pcaX <- dudi.pca(data, center = center, scale = scale, scannf = FALSE, nf=maxRank)
		  pcs <- pcaX$li
		  # then perform GUIDED random forest clustering
		  rf.fit <- randomForest(x = pcs, y = model, ntree = 10000, proximity = TRUE, oob.prox = TRUE)
		  hclust.rf <- hclust(as.dist(1-rf.fit$proximity), method = "ward.D2")
		  rf.cluster = cutree(hclust.rf, k=k)
		  dapc.rf <- dapc(data, rf.cluster, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.var, n.da=10000, center=center, scale=scale)


		  # write RF group assignment to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.supervised.grp.txt", sep=""), x = rf.cluster,
					  row.names = FALSE, col.names = FALSE)

		  # write DAPC variable contributions to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.supervised.var.contr.txt", sep=""), x = dapc.rf$var.contr)

		  # write DAPC variable loadings to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.supervised.var.load.txt", sep=""), x = dapc.rf$var.load)

		  # write DAPC individual coordinates to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.supervised.ind.coord.txt", sep=""), x = dapc.rf$ind.coord,
					  row.names=FALSE)

		  # write DAPC assignment probabilities to file
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, "_RF.supervised.posteriors.txt", sep=""), x = dapc.rf$posterior,
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
		  d2 <- data.frame(1:nrow(dapc.rf$posterior),
						   1:nrow(dapc.rf$posterior),
						   rep("(x)", nrow(dapc.rf$posterior)),
						   model,
						   rep(":", nrow(dapc.rf$posterior)),
						   dapc.rf$posterior)
		  colnames(d2) <- colnames(d1)
		  d <- rbind(d1, d2)
		  # create CLUMPP indfile
		  write.table(file = paste("m", modelNumber, "_perVar-", perc.var, ".RF.supervised.indfile", sep=""), x=d,
					  quote=FALSE, col.names = FALSE, row.names=FALSE)

		  # create CLUMPP paramfile
		  f <- paste("m", modelNumber, "_perVar-", perc.var, ".RF.supervised.paramfile", sep="")
		  file <- file(f, "wb")
		  cat("DATATYPE 0", file=file, sep="\n")
		  cat(paste("INDFILE m", modelNumber, "_perVar-", perc.var, ".RF.supervised.indfile", sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("OUTFILE m", modelNumber, "_perVar-", perc.var, ".RF.supervised.outfile", sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("MISCFILE m", modelNumber, "_perVar-", perc.var, ".RF.supervised.miscfile", sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("K ", ncol(pop.assign), sep=""), file=file, append=TRUE, sep="\n")
		  cat(paste("C ", nrow(pop.assign), sep=""), file=file, append=TRUE, sep="\n")
		  cat("R 2\nM 1\nW 0\nS 2\nPRINT_PERMUTED_DATA 0\nPRINT_EVERY_PERM 0\nPRINT_RANDOM_INPUTORDER 0\nOVERRIDE_WARNINGS 0\nORDER_BY_RUN 1", file=file, append=TRUE, sep="\n")

		  shell(paste(CLUMPP_exe, f, sep=" "))

		  close(file)
	  }
  }
}


