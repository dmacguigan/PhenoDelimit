# create scatter plot or density plot of LD axes
plot_discriminant_axes <- function(wd, clumpp.wd,
                                   sample.plot.groups = NULL, sample.plot.groups.order = NULL,
                                   best.perc.var, best.model.number,
                                   plot.type, plot.width, plot.height, colors, shapes,
                                   x.axis=NULL, y.axis=NULL){

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
      svglite(file = paste("DA", x.axis, "-DA", y.axis, "_scatter", ".svg", sep=""), width = plot.width, height = plot.height)
      print(plot)
      dev.off()
    } else if(plot.type =="pdf"){
      pdf(file = paste("DA", x.axis, "-DA", y.axis, "_scatter", ".pdf", sep=""), width = plot.width, height = plot.height)
      print(plot)
      dev.off()
    } else {
      png(file = paste("DA", x.axis, "-DA", y.axis, "_scatter", ".png", sep=""), units="in", res=300,width = plot.width, height = plot.height)
      print(plot)
      dev.off()
    }

    return(plot)

    return(plot)
  } else {
    p <- data.frame(ind.coord)
    p$cluster <- factor(kmeans.grp)
    p$group <- factor(sample.plot.groups, levels=sample.plot.groups.order)

    plot <- ggplot(p, aes(x=LD1, color=group, fill=group)) +
      geom_density(alpha=0.5) +
      scale_color_manual(values = colors, labels = sample.plot.groups.order, name = "Group") +
      scale_fill_manual(values = colors, guide="none") +
      theme_minimal() +
      xlab("Discriminant axis 1")

    if(plot.type =="svg"){
      svglite(file = "DA1_density.svg", width = plot.width, height = plot.height)
      print(plot)
      dev.off()
    } else if(plot.type =="pdf"){
      pdf(file = "DA1_density.pdf", width = plot.width, height = plot.height)
      print(plot)
      dev.off()
    } else {
      png(file = "DA1_density.png", units="in", res=300,width = plot.width, height = plot.height)
      print(plot)
      dev.off()
    }

    return(plot)
  }

}




