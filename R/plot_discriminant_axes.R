# create scatter plot or density plot of LD axes
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
		hulls <- ddply(p, "group", find_hull)

		p$cluster <- factor(kmeans.grp[,1])

		plot <- ggplot() +
		  geom_polygon(data=hulls, aes_string(x=paste("LD", x.axis, sep=""), y=paste("LD", y.axis, sep=""),
											  color="group", fill="group"), alpha=0.2) +
		  geom_point(data=p, aes_string(x=paste("LD", x.axis, sep=""), y=paste("LD", y.axis, sep=""),
										color="group", fill="group", shape="cluster"), cex=3, alpha=0.8) +
		  scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  scale_fill_manual(values = colors, guide="none") +
		  scale_shape_manual(values = shapes, labels = paste(1:(ncol(ind.coord)+1)), name = "K-means\nCluster") +
		  theme(panel.background = element_rect(fill = 'white', colour = 'black'), legend.box = "horizontal") +
		  theme_minimal() +
		  xlab(paste("Discriminant axis ", x.axis)) +
		  ylab(paste("Discriminant axis ", y.axis))

		if(plot.type =="svg"){
		  svglite(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_scatter", ".svg", sep=""), width = plot.width, height = plot.height)
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

		plot <- ggplot(p, aes(x=LD1, color=group, fill=group)) +
		  geom_density(alpha=0.5) +
		  scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  scale_fill_manual(values = colors, guide="none") +
		  theme_minimal() +
		  xlab("Discriminant axis 1")

		if(plot.type =="svg"){
		  svglite(file = paste0("m", best.model.number, "_DA1_density.svg"), width = plot.width, height = plot.height)
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
		hulls <- ddply(p, "group", find_hull)

		plot <- ggplot() +
		  geom_polygon(data=hulls, aes_string(x=paste("LD", x.axis, sep=""), y=paste("LD", y.axis, sep=""),
											  color="group", fill="group"), alpha=0.2) +
		  geom_point(data=p, aes_string(x=paste("LD", x.axis, sep=""), y=paste("LD", y.axis, sep=""),
										color="group", fill="group"), cex=3, alpha=0.8) +
		  scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  scale_fill_manual(values = colors, guide="none") +
		  theme(panel.background = element_rect(fill = 'white', colour = 'black'), legend.box = "horizontal") +
		  theme_minimal() +
		  xlab(paste("Discriminant axis ", x.axis)) +
		  ylab(paste("Discriminant axis ", y.axis))

		if(plot.type =="svg"){
		  svglite(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_scatter", ".apriori.svg", sep=""), width = plot.width, height = plot.height)
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

		plot <- ggplot(p, aes(x=LD1, color=group, fill=group)) +
		  geom_density(alpha=0.5) +
		  scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  scale_fill_manual(values = colors, guide="none") +
		  theme_minimal() +
		  xlab("Discriminant axis 1")

		if(plot.type =="svg"){
		  svglite(file = paste0("m", best.model.number, "_DA1_density.apriori.svg"), width = plot.width, height = plot.height)
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
		hulls <- ddply(p, "group", find_hull)

		p$cluster <- factor(rf.grp[,1])

		plot <- ggplot() +
		  geom_polygon(data=hulls, aes_string(x=paste("LD", x.axis, sep=""), y=paste("LD", y.axis, sep=""),
											  color="group", fill="group"), alpha=0.2) +
		  geom_point(data=p, aes_string(x=paste("LD", x.axis, sep=""), y=paste("LD", y.axis, sep=""),
										color="group", fill="group", shape="cluster"), cex=3, alpha=0.8) +
		  scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  scale_fill_manual(values = colors, guide="none") +
		  scale_shape_manual(values = shapes, labels = paste(1:(ncol(ind.coord)+1)), name = "RF\nCluster") +
		  theme(panel.background = element_rect(fill = 'white', colour = 'black'), legend.box = "horizontal") +
		  theme_minimal() +
		  xlab(paste("Discriminant axis ", x.axis)) +
		  ylab(paste("Discriminant axis ", y.axis))

		if(plot.type =="svg"){
		  svglite(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_RF.scatter", ".svg", sep=""), width = plot.width, height = plot.height)
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

		plot <- ggplot(p, aes(x=LD1, color=group, fill=group)) +
		  geom_density(alpha=0.5) +
		  scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  scale_fill_manual(values = colors, guide="none") +
		  theme_minimal() +
		  xlab("Discriminant axis 1")

		if(plot.type =="svg"){
		  svglite(file = paste0("m", best.model.number, "_DA1_RF.density.svg"), width = plot.width, height = plot.height)
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
		hulls <- ddply(p, "group", find_hull)

		p$cluster <- factor(rf.grp[,1])

		plot <- ggplot() +
		  geom_polygon(data=hulls, aes_string(x=paste("LD", x.axis, sep=""), y=paste("LD", y.axis, sep=""),
											  color="group", fill="group"), alpha=0.2) +
		  geom_point(data=p, aes_string(x=paste("LD", x.axis, sep=""), y=paste("LD", y.axis, sep=""),
										color="group", fill="group", shape="cluster"), cex=3, alpha=0.8) +
		  scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  scale_fill_manual(values = colors, guide="none") +
		  scale_shape_manual(values = shapes, labels = paste(1:(ncol(ind.coord)+1)), name = "RF\nCluster") +
		  theme(panel.background = element_rect(fill = 'white', colour = 'black'), legend.box = "horizontal") +
		  theme_minimal() +
		  xlab(paste("Discriminant axis ", x.axis)) +
		  ylab(paste("Discriminant axis ", y.axis))

		if(plot.type =="svg"){
		  svglite(file = paste("m", best.model.number, "_DA", x.axis, "-DA", y.axis, "_RF.supervised.scatter", ".svg", sep=""), width = plot.width, height = plot.height)
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

		plot <- ggplot(p, aes(x=LD1, color=group, fill=group)) +
		  geom_density(alpha=0.5) +
		  scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
		  scale_fill_manual(values = colors, guide="none") +
		  theme_minimal() +
		  xlab("Discriminant axis 1")

		if(plot.type =="svg"){
		  svglite(file = paste0("m", best.model.number, "_DA1_RF.supervised.density.svg"), width = plot.width, height = plot.height)
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




