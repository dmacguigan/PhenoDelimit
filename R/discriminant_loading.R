#' DAPC loading plots
#' 
#' Make loading plot of variables on selected discriminant axis. Saves plot as file (svg or pdf).
#' Also returns file with variable contributions and loadings for that axis.
#' 
#' @param wd working directory to create plots
#' @param clumpp.wd directory containing CLUMPP results
#' @param best.perc.var which percent retained variance to plot? numeric
#' @param best.model.number which model number to plot? numeric
#' @param plot.type save plot as "pdf", "svg", or "png"
#' @param plot.width width of plot in inches
#' @param plot.height height of plot in inches
#' @param axis which discriminant axis to examine? numeric
#' @param apriori do you wish to plot results from apriori individual assignment to species/populations/clusters, or results from assignment using using k-means clustering? TRUE or FALSE
#' 
#' @export

discriminant_loading <- function(wd, clumpp.wd,
                                 best.perc.var, best.model.number,
                                 plot.type, plot.width, plot.height,
                                 axis, clust.method="kmeans", apriori=FALSE){

  if(apriori == FALSE && clust.method == "kmeans"){
	  
	  # read DAPC variable contributions from file
	  var.contr <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_var.contr.txt", sep=""),
							  header = TRUE)

	  # read DAPC variable loadings from file
	  var.load <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_var.load.txt", sep=""),
							 header=TRUE)

	  setwd(wd)

	  if(plot.type =="svg"){
		svglite::svglite(file = paste("m", best.model.number, "_DA_", axis, "_loading.svg", sep=""), width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  } else if(plot.type =="pdf"){
		grDevices::pdf(file = paste("m", best.model.number, "_DA_", axis, "_loading.pdf", sep=""), width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  } else {
		grDevices::png(file = paste("m", best.model.number, "_DA_", axis, "_loading.png", sep=""), units="in", res=300, width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  }
	   
	  adegenet::loadingplot(var.contr, axis=axis, threshold=0)

	  d <- var.contr[,axis]
	  d2 <- var.load[,axis]
	  t <- cbind(d,d2)
	  t <- as.data.frame(t)
	  row.names(t) <- row.names(var.contr)
	  t <- t[order(t$d, decreasing = TRUE),]
	  colnames(t) <- c("contribution", "loading")
	  write.table(x = t, file = paste("m", best.model.number, "_DA_", axis, "_loading.txt", sep=""), quote=FALSE, row.names = TRUE)
  }
  
  else if(apriori == TRUE && clust.method == "kmeans"){
  	# read DAPC variable contributions from file
	  var.contr <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_var.contr.apriori.txt", sep=""),
							  header = TRUE)

	  # read DAPC variable loadings from file
	  var.load <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_var.load.apriori.txt", sep=""),
							 header=TRUE)

	  setwd(wd)

	  if(plot.type =="svg"){
		svglite::svglite(file = paste("m", best.model.number, "_DA_", axis, "_loading.apriori.svg", sep=""), width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  } else if(plot.type =="pdf"){
		grDevices::pdf(file = paste("m", best.model.number, "_DA_", axis, "_loading.apriori.pdf", sep=""), width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  } else {
		grDevices::png(file = paste("m", best.model.number, "_DA_", axis, "_loading.apriori.png", sep=""), units="in", res=300, width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  }
	   
	  adegenet::loadingplot(var.contr, axis=axis, threshold=0)

	  d <- var.contr[,axis]
	  d2 <- var.load[,axis]
	  t <- cbind(d,d2)
	  t <- as.data.frame(t)
	  row.names(t) <- row.names(var.contr)
	  t <- t[order(t$d, decreasing = TRUE),]
	  colnames(t) <- c("contribution", "loading")
	  write.table(x = t, file = paste("m", best.model.number, "_DA_", axis, "_loading.apriori.txt", sep=""), quote=FALSE, row.names = TRUE)
  
  }
  
  else if(apriori == FALSE && clust.method == "randomforest"){
	  
	  # read DAPC variable contributions from file
	  var.contr <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_RF.var.contr.txt", sep=""),
							  header = TRUE)

	  # read DAPC variable loadings from file
	  var.load <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_RF.var.load.txt", sep=""),
							 header=TRUE)

	  setwd(wd)

	  if(plot.type =="svg"){
		svglite::svglite(file = paste("m", best.model.number, "_DA_", axis, "_RF.loading.svg", sep=""), width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  } else if(plot.type =="pdf"){
		grDevices::pdf(file = paste("m", best.model.number, "_DA_", axis, "_RF.loading.pdf", sep=""), width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  } else {
		grDevices::png(file = paste("m", best.model.number, "_DA_", axis, "_RF.loading.png", sep=""), units="in", res=300, width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  }
	   
	  adegenet::loadingplot(var.contr, axis=axis, threshold=0)

	  d <- var.contr[,axis]
	  d2 <- var.load[,axis]
	  t <- cbind(d,d2)
	  t <- as.data.frame(t)
	  row.names(t) <- row.names(var.contr)
	  t <- t[order(t$d, decreasing = TRUE),]
	  colnames(t) <- c("contribution", "loading")
	  write.table(x = t, file = paste("m", best.model.number, "_DA_", axis, "_RF.loading.txt", sep=""), quote=FALSE, row.names = TRUE)
  }
  
    else if(apriori == TRUE && clust.method == "randomforest"){
	  
	  # read DAPC variable contributions from file
	  var.contr <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_RF.supervised.var.contr.txt", sep=""),
							  header = TRUE)

	  # read DAPC variable loadings from file
	  var.load <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_RF.supervised.var.load.txt", sep=""),
							 header=TRUE)

	  setwd(wd)

	  if(plot.type =="svg"){
		svglite::svglite(file = paste("m", best.model.number, "_DA_", axis, "_RF.supervised.loading.svg", sep=""), width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  } else if(plot.type =="pdf"){
		grDevices::pdf(file = paste("m", best.model.number, "_DA_", axis, "_RF.supervised.loading.pdf", sep=""), width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  } else {
		grDevices::png(file = paste("m", best.model.number, "_DA_", axis, "_RF.supervised.loading.png", sep=""), units="in", res=300, width = plot.width, height = plot.height)
		adegenet::loadingplot(var.contr, axis=axis, threshold=0)
		dev.off()
	  }
	   
	  adegenet::loadingplot(var.contr, axis=axis, threshold=0)

	  d <- var.contr[,axis]
	  d2 <- var.load[,axis]
	  t <- cbind(d,d2)
	  t <- as.data.frame(t)
	  row.names(t) <- row.names(var.contr)
	  t <- t[order(t$d, decreasing = TRUE),]
	  colnames(t) <- c("contribution", "loading")
	  write.table(x = t, file = paste("m", best.model.number, "_DA_", axis, "_RF.supervised.loading.txt", sep=""), quote=FALSE, row.names = TRUE)
  }
}


