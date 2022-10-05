# function to create CLUMPP indfiles for model assignment vs K-means clusters DAPC probs
# also allows the user to run DAPC with a priori species/group assignment
dapc_clumpp <- function(wd, data, n.groups, model.numbers, models, perc.var, scale=TRUE, center=TRUE, apriori=FALSE){
  setwd(wd)
  for(i in 1:length(n.groups)){
    for(j in 1:length(perc.var)){
      run_dapc_clumpp(data=data, modelNumber=model.numbers[i], k=n.groups[i], model=as.factor(models[,i]), perc.var=perc.var[j], center=center, scale=scale, apriori=apriori)
    }
  }
}

run_dapc_clumpp <- function(data, modelNumber, k, model, perc.var, scale, center, apriori){
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

	  shell(paste("CLUMPP", f, sep=" "))
	  
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
}


