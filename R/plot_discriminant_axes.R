#' Make scatter or density plots of discriminant axes
#' 
#' Make density plot (2 clusters) or scatter plot (>2 clusters) of discriminant axes, as specified by the user. Returns plot and saves plot as file (svg or pdf).
#' 
#' @param wd directory to create plots
#' @param clumpp.wd directory containing CLUMPP results
#' @param sample.names vector of sample names, same order as original data, numeric or character
#' @param sample.plot.groups vector of sample groups, same order as original data, numeric or character
#' @param sample.plot.groups.order vector of sample groups in desired order, numeric or character
#' @param best.perc.var which percent retained variance to plot? numeric
#' @param best.model.number which model number to plot? numeric
#' @param plot.type save plot as "pdf", "svg", or "png"
#' @param plot.width width of plot in inches
#' @param plot.height height of plot in inches
#' @param colors vector of colors for groups
#' @param shapes vector of shapes for K-means clusters
#' @param x.axis which discriminant axis to show on x axis? If only two clusters, discriminant axis 1 is automatically used. numeric
#' @param y.axis which discriminant axis to show on y axis? If only two clusters, this argument is not used. numeric
#' @param apriori do you wish to plot results from apriori individual assignment to species/populations/clusters, or results from assignment using using k-means clustering? TRUE or FALSE
#' 
#' @export
#' 
plot_discriminant_axes <- function(wd, clumpp.wd,
                                   sample.plot.groups = NULL, sample.plot.groups.order = NULL,
                                   best.perc.var, best.model.number,
                                   plot.type, plot.width, plot.height, colors, shapes,
                                   x.axis=NULL, y.axis=NULL,
                                   apriori=FALSE, clust.method="kmeans"){

  if(apriori == FALSE && clust.method == "kmeans"){
	  # read DAPC individual coordinates from file
	  ind.coord <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_ind.coord.txt", sep=""),
							  header=TRUE)

	  # read K-means group assignment from file
	  kmeans.grp <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_Kmeans.grp.txt", sep=""),
							   header=FALSE)

	  setwd(wd)

	  if(ncol(ind.coord) > 1){
		p <- data.frame(ind.coord)
		p$group <- factor(sample.plot.groups, levels=sample.plot.groups.order)

		find_hull <- function(df) df[chull(df[,x.axis], df[,y.axis]), ]
		hulls <- plyr::ddply(p, "group", find_hull)

		p$cluster <- factor(kmeans.grp[,1])

		plot <- ggplot2::ggplot() +
		  ggplot2::geom_polygon(data = hulls, 
                     ggplot2::aes(x = .data[[paste("LD", x.axis, sep="")]], 
                                 y = .data[[paste("LD", y.axis, sep="")]], 
                                 color = .data[["group"]], 
                                 fill = .data[["group"]]), 
                     alpha = 0.2) +
			ggplot2::geom_point(data = p, 
                   ggplot2::aes(x = .data[[paste("LD", x.axis, sep="")]], 
                               y = .data[[paste("LD", y.axis, sep="")]], 
                               color = .data[["group"]], 
                               fill = .data[["group"]], 
                               shape = .data[["cluster"]]), 
                   size = 3, alpha = 0.8) +
		  ggplot2::scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  ggplot2::scale_fill_manual(values = colors, guide="none") +
		  ggplot2::scale_shape_manual(values = shapes, labels = paste(1:(ncol(ind.coord)+1)), name = "K-means\nCluster") +
		  ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'black'), legend.box = "horizontal") +
		  ggplot2::theme_minimal() +
		  ggplot2::xlab(paste("Discriminant axis ", x.axis)) +
		  ggplot2::ylab(paste("Discriminant axis ", y.axis))

		if(plot.type =="svg"){
		  svglite::svglite(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_scatter", ".svg", sep=""), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else if(plot.type =="pdf"){
		  pdf(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_scatter", ".pdf", sep=""), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else {
		  png(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_scatter", ".png", sep=""), units="in", res=300,width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		}

		return(plot)

		return(plot)
	  } else {
		p <- data.frame(ind.coord)
		p$cluster <- factor(kmeans.grp[,1])
		p$group <- factor(sample.plot.groups, levels=sample.plot.groups.order)

		plot <- ggplot2::ggplot(p, ggplot2::aes(x=LD1, color=group, fill=group)) +
		  ggplot2::geom_density(alpha=0.5) +
		  ggplot2::scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  ggplot2::scale_fill_manual(values = colors, guide="none") +
		  ggplot2::theme_minimal() +
		  ggplot2::xlab("Discriminant axis 1")

		if(plot.type =="svg"){
		  svglite::svglite(file = paste0("m", best.model.number, "_DA1_density.svg"), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else if(plot.type =="pdf"){
		  pdf(file = paste0("m", best.model.number, "_DA1_density.pdf"), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else {
		  png(file = paste0("m", best.model.number, "_DA1_density.png"), units="in", res=300,width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		}

		return(plot)
	  }
  }else if(apriori == TRUE && clust.method == "kmeans") {
  	  # read DAPC individual coordinates from file
	  ind.coord <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_ind.coord.apriori.txt", sep=""),
							  header=TRUE)


	  setwd(wd)

	  if(ncol(ind.coord) > 1){
		p <- data.frame(ind.coord)
		p$group <- factor(sample.plot.groups, levels=sample.plot.groups.order)

		find_hull <- function(df) df[chull(df[,x.axis], df[,y.axis]), ]
		hulls <- plyr::ddply(p, "group", find_hull)

		plot <- ggplot2::ggplot() +
		  ggplot2::geom_polygon(data = hulls, 
                     ggplot2::aes(x = .data[[paste("LD", x.axis, sep="")]], 
                                 y = .data[[paste("LD", y.axis, sep="")]], 
                                 color = .data[["group"]], 
                                 fill = .data[["group"]]), 
                     alpha = 0.2) +
			ggplot2::geom_point(data = p, 
                   ggplot2::aes(x = .data[[paste("LD", x.axis, sep="")]], 
                               y = .data[[paste("LD", y.axis, sep="")]], 
                               color = .data[["group"]], 
                               fill = .data[["group"]]), 
                   size = 3, alpha = 0.8) +
		  ggplot2::scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  ggplot2::scale_fill_manual(values = colors, guide="none") +
		  ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'black'), legend.box = "horizontal") +
		  ggplot2::theme_minimal() +
		  ggplot2::xlab(paste("Discriminant axis ", x.axis)) +
		  ggplot2::ylab(paste("Discriminant axis ", y.axis))

		if(plot.type =="svg"){
		  svglite::svglite(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_scatter", ".apriori.svg", sep=""), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else if(plot.type =="pdf"){
		  pdf(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_scatter", ".apriori.pdf", sep=""), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else {
		  png(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_scatter", ".apriori.png", sep=""), units="in", res=300,width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		}

		return(plot)

		return(plot)
	  } else {
		p <- data.frame(ind.coord)
		p$group <- factor(sample.plot.groups, levels=sample.plot.groups.order)

		plot <- ggplot2::ggplot(p, ggplot2::aes(x=LD1, color=group, fill=group)) +
		  ggplot2::geom_density(alpha=0.5) +
		  ggplot2::scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  ggplot2::scale_fill_manual(values = colors, guide="none") +
		  ggplot2::theme_minimal() +
		  ggplot2::xlab("Discriminant axis 1")

		if(plot.type =="svg"){
		  svglite::svglite(file = paste0("m", best.model.number, "_DA1_density.apriori.svg"), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else if(plot.type =="pdf"){
		  pdf(file = paste0("m", best.model.number, "_DA1_density.apriori.pdf"), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else {
		  png(file = paste0("m", best.model.number, "_DA1_density.apriori.png"), units="in", res=300,width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		}

		return(plot)
	  }
  }else if(apriori == FALSE && clust.method == "randomforest"){
	  # read DAPC individual coordinates from file
	  ind.coord <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_RF.ind.coord.txt", sep=""),
							  header=TRUE)

	  # read RF group assignment from file
	  rf.grp <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_RF.grp.txt", sep=""),
							   header=FALSE)

	  setwd(wd)

	  if(ncol(ind.coord) > 1){
		p <- data.frame(ind.coord)
		p$group <- factor(sample.plot.groups, levels=sample.plot.groups.order)

		find_hull <- function(df) df[chull(df[,x.axis], df[,y.axis]), ]
		hulls <- plyr::ddply(p, "group", find_hull)

		p$cluster <- factor(rf.grp[,1])

		plot <- ggplot2::ggplot() +
			ggplot2::geom_polygon(data = hulls, 
                     ggplot2::aes(x = .data[[paste("LD", x.axis, sep="")]], 
                                 y = .data[[paste("LD", y.axis, sep="")]], 
                                 color = .data[["group"]], 
                                 fill = .data[["group"]]), 
                     alpha = 0.2) +
			ggplot2::geom_point(data = p, 
                   ggplot2::aes(x = .data[[paste("LD", x.axis, sep="")]], 
                               y = .data[[paste("LD", y.axis, sep="")]], 
                               color = .data[["group"]], 
                               fill = .data[["group"]], 
                               shape = .data[["cluster"]]), 
                   size = 3, alpha = 0.8) +
		  ggplot2::scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  ggplot2::scale_fill_manual(values = colors, guide="none") +
		  ggplot2::scale_shape_manual(values = shapes, labels = paste(1:(ncol(ind.coord)+1)), name = "RF\nCluster") +
		  ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'black'), legend.box = "horizontal") +
		  ggplot2::theme_minimal() +
		  ggplot2::xlab(paste("Discriminant axis ", x.axis)) +
		  ggplot2::ylab(paste("Discriminant axis ", y.axis))

		if(plot.type =="svg"){
		  svglite::svglite(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_RF.scatter", ".svg", sep=""), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else if(plot.type =="pdf"){
		  pdf(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_RF.scatter", ".pdf", sep=""), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else {
		  png(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_RF.scatter", ".png", sep=""), units="in", res=300,width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		}

		return(plot)

		return(plot)
	  } else {
		p <- data.frame(ind.coord)
		p$cluster <- factor(rf.grp[,1])
		p$group <- factor(sample.plot.groups, levels=sample.plot.groups.order)

		plot <- ggplot2::ggplot(p, ggplot2::aes(x=LD1, color=group, fill=group)) +
		  ggplot2::eom_density(alpha=0.5) +
		  ggplot2::scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  ggplot2::scale_fill_manual(values = colors, guide="none") +
		  ggplot2::theme_minimal() +
		  ggplot2::xlab("Discriminant axis 1")

		if(plot.type =="svg"){
		  svglite::svglite(file = paste0("m", best.model.number, "_DA1_RF.density.svg"), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else if(plot.type =="pdf"){
		  pdf(file = paste0("m", best.model.number, "_DA1_RF.density.pdf"), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else {
		  png(file = paste0("m", best.model.number, "_DA1_RF.density.png"), units="in", res=300,width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		}

		return(plot)
	  }
  }else if(apriori == TRUE && clust.method == "randomforest") {
	  # read DAPC individual coordinates from file
	  ind.coord <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_RF.supervised.ind.coord.txt", sep=""),
							  header=TRUE)

	  # read RF group assignment from file
	  rf.grp <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_RF.supervised.grp.txt", sep=""),
							   header=FALSE)

	  setwd(wd)

	  if(ncol(ind.coord) > 1){
		p <- data.frame(ind.coord)
		p$group <- factor(sample.plot.groups, levels=sample.plot.groups.order)

		find_hull <- function(df) df[chull(df[,x.axis], df[,y.axis]), ]
		hulls <- plyr::ddply(p, "group", find_hull)

		p$cluster <- factor(rf.grp[,1])

		plot <- ggplot2::ggplot() +
			ggplot2::geom_polygon(data = hulls, 
                     ggplot2::aes(x = .data[[paste("LD", x.axis, sep="")]], 
                                 y = .data[[paste("LD", y.axis, sep="")]], 
                                 color = .data[["group"]], 
                                 fill = .data[["group"]]), 
                     alpha = 0.2) +
			ggplot2::geom_point(data = p, 
                   ggplot2::aes(x = .data[[paste("LD", x.axis, sep="")]], 
                               y = .data[[paste("LD", y.axis, sep="")]], 
                               color = .data[["group"]], 
                               fill = .data[["group"]]), 
                   size = 3, alpha = 0.8) +
		  ggplot2::scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  ggplot2::scale_fill_manual(values = colors, guide="none") +
		  ggplot2::scale_shape_manual(values = shapes, labels = paste(1:(ncol(ind.coord)+1)), name = "RF\nCluster") +
		  ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'black'), legend.box = "horizontal") +
		  ggplot2::theme_minimal() +
		  ggplot2::xlab(paste("Discriminant axis ", x.axis)) +
		  ggplot2::ylab(paste("Discriminant axis ", y.axis))

		if(plot.type =="svg"){
		  svglite::svglite(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_RF.supervised.scatter", ".svg", sep=""), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else if(plot.type =="pdf"){
		  pdf(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_RF.supervised.scatter", ".pdf", sep=""), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else {
		  png(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_RF.supervised.scatter", ".png", sep=""), units="in", res=300,width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		}

		return(plot)

		return(plot)
	  } else {
		p <- data.frame(ind.coord)
		p$cluster <- factor(rf.grp[,1])
		p$group <- factor(sample.plot.groups, levels=sample.plot.groups.order)

		plot <- ggplot2::ggplot(p, ggplot2::aes(x=LD1, color=group, fill=group)) +
		  ggplot2::geom_density(alpha=0.5) +
		  ggplot2::scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  ggplot2::scale_fill_manual(values = colors, guide="none") +
		  ggplot2::theme_minimal() +
		  ggplot2::xlab("Discriminant axis 1")

		if(plot.type =="svg"){
		  svglite::svglite(file = paste0("m", best.model.number, "_DA1_RF.supervised.density.svg"), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else if(plot.type =="pdf"){
		  pdf(file = paste0("m", best.model.number, "_DA1_RF.supervised.density.pdf"), width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		} else {
		  png(file = paste0("m", best.model.number, "_DA1_RF.supervised.density.png"), units="in", res=300,width = plot.width, height = plot.height)
		  print(plot)
		  dev.off()
		}

		return(plot)
	  }
	}

}




