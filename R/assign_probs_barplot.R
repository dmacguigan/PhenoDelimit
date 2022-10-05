# make barplot of best result posterior probabilities
assign_probs_barplot <- function(wd, clumpp.wd, sample.names,
                                sample.plot.groups = NULL, sample.plot.groups.order = NULL,
                                sample.order = NULL,
                                best.perc.var, best.model.number,
                                plot.type, plot.width, plot.height, plot.name, colors, border.color,
								apriori=FALSE){
  
  if(apriori == FALSE){
	  # read in posteriors from clumpp_prep function
	  dapc.data <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_posteriors.txt", sep=""))

	  # order the data for plotting
	  if(!(is.null(sample.order))){
		# create order vector based on sample.plot.groups.order
		pop.order <- factor(sample.plot.groups, levels=sample.plot.groups.order)
		plotData <- readData_sampleOrder(dapc.data = dapc.data, sample.names = sample.names, pop.order = pop.order, sample.order = sample.order)

	  } else if(!(is.null(sample.plot.groups))){
		# create order vector based on sample.plot.groups.order
		pop.order <- factor(sample.plot.groups, levels=sample.plot.groups.order)
		plotData <- readData(dapc.data = dapc.data,  sample.names = sample.names, pop.order = pop.order)

	  }

	  # get positions of breaks between populations
	  pop.order.ordered <- sort(pop.order)
	  pop.lines <- cumsum(table(pop.order.ordered))
	  pop.labels <- (pop.lines - table(pop.order.ordered)) + (table(pop.order.ordered)/2)

	  # make barplots
	  setwd(wd)
	  if(plot.type == "svg"){
		svglite(file=paste(plot.name, ".svg", sep=""), width=plot.width, height=plot.height)
		  par(mgp=c(1.5,1,0), mar = c(4.1, 3.1, 1, 1))
		  barplot(plotData, col=colors, border=border.color,
				  space=0, axes=F, axisname=FALSE, las=2, cex.names=0.75,
				  ylab  = "Assignment Probability", xlab = "Group")
		  for(i in 1:(length(pop.lines)-1)){
			abline(v=pop.lines[i], lwd=2, lty=2, col="black") # set at divisions between populations
		  }
		  axis(2, las=2, cex=0.75, pos=0)
		  axis(1, at=c(seq(from = 0.5, to = ncol(plotData + 0.5), by = 1)), labels=sample.names,
			   tick=FALSE, las=2, cex.axis=0.3, adj=0, pos=0, line=0, mgp = c(3, 0.1, 0))
		  axis(1, at=pop.labels, labels=names(table(pop.order.ordered)),
			   tick=FALSE, las=0, cex.axis=0.8, adj=0, pos=0)
		dev.off()
	  } else if(plot.type == "pdf"){
		pdf(file=paste(plot.name, ".pdf", sep=""), width=plot.width, height=plot.height)
		  par(mgp=c(1.5,1,0), mar = c(4.1, 3.1, 1, 1))
		  barplot(plotData, col=colors, border=border.color,
				  space=0, axes=F, axisname=FALSE, las=2, cex.names=0.75,
				  ylab  = "Assignment Probability", xlab = "Group")
		  for(i in 1:(length(pop.lines)-1)){
			abline(v=pop.lines[i], lwd=2, lty=2, col="black") # set at divisions between populations
		  }
		  axis(2, las=2, cex=0.75, pos=0)
		  axis(1, at=c(seq(from = 0.5, to = ncol(plotData + 0.5), by = 1)), labels=sample.names,
			   tick=FALSE, las=2, cex.axis=0.3, adj=0, pos=0, line=0, mgp = c(3, 0.1, 0))
		  axis(1, at=pop.labels, labels=names(table(pop.order.ordered)),
			   tick=FALSE, las=0, cex.axis=0.8, adj=0, pos=0)
		dev.off()
	  } else {
		png(file=paste(plot.name, ".png", sep=""), units="in", res=300, width=plot.width, height=plot.height)
		par(mgp=c(1.5,1,0), mar = c(4.1, 3.1, 1, 1))
		barplot(plotData, col=colors, border=border.color,
				space=0, axes=F, axisname=FALSE, las=2, cex.names=0.75,
				ylab  = "Assignment Probability", xlab = "Group")
		for(i in 1:(length(pop.lines)-1)){
		  abline(v=pop.lines[i], lwd=2, lty=2, col="black") # set at divisions between populations
		}
		axis(2, las=2, cex=0.75, pos=0)
		axis(1, at=c(seq(from = 0.5, to = ncol(plotData + 0.5), by = 1)), labels=sample.names,
			 tick=FALSE, las=2, cex.axis=0.3, adj=0, pos=0, line=0, mgp = c(3, 0.1, 0))
		axis(1, at=pop.labels, labels=names(table(pop.order.ordered)),
			 tick=FALSE, las=0, cex.axis=0.8, adj=0, pos=0)
		dev.off()
	  }

	  par(mgp=c(1.5,1,0), mar = c(4.1, 3.1, 1, 1))
	  barplot(plotData, col=colors, border=border.color,
			  space=0, axes=F, axisname=FALSE, las=2, cex.names=0.75,
			  ylab  = "Assignment Probability", xlab = "Group")
	  for(i in 1:(length(pop.lines)-1)){
		abline(v=pop.lines[i], lwd=2, lty=2, col="black") # set at divisions between populations
	  }
	  axis(2, las=2, cex=0.75, pos=0)
	  axis(1, at=c(seq(from = 0.5, to = ncol(plotData + 0.5), by = 1)), labels=sample.names,
		   tick=FALSE, las=2, cex.axis=0.3, adj=0, pos=0, line=0, mgp = c(3, 0.1, 0))
	  axis(1, at=pop.labels, labels=names(table(pop.order.ordered)),
		   tick=FALSE, las=0, cex.axis=0.8, adj=0, pos=0)
	  title()
  }
  else{
  	  # read in posteriors from clumpp_prep function
	  dapc.data <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_posteriors.apriori.txt", sep=""))

	  # order the data for plotting
	  if(!(is.null(sample.order))){
		# create order vector based on sample.plot.groups.order
		pop.order <- factor(sample.plot.groups, levels=sample.plot.groups.order)
		plotData <- readData_sampleOrder(dapc.data = dapc.data, sample.names = sample.names, pop.order = pop.order, sample.order = sample.order)

	  } else if(!(is.null(sample.plot.groups))){
		# create order vector based on sample.plot.groups.order
		pop.order <- factor(sample.plot.groups, levels=sample.plot.groups.order)
		plotData <- readData(dapc.data = dapc.data,  sample.names = sample.names, pop.order = pop.order)

	  }

	  # get positions of breaks between populations
	  pop.order.ordered <- sort(pop.order)
	  pop.lines <- cumsum(table(pop.order.ordered))
	  pop.labels <- (pop.lines - table(pop.order.ordered)) + (table(pop.order.ordered)/2)

	  # make barplots
	  setwd(wd)
	  if(plot.type == "svg"){
		svglite(file=paste(plot.name, ".apriori.svg", sep=""), width=plot.width, height=plot.height)
		  par(mgp=c(1.5,1,0), mar = c(4.1, 3.1, 1, 1))
		  barplot(plotData, col=colors, border=border.color,
				  space=0, axes=F, axisname=FALSE, las=2, cex.names=0.75,
				  ylab  = "Assignment Probability", xlab = "Group")
		  for(i in 1:(length(pop.lines)-1)){
			abline(v=pop.lines[i], lwd=2, lty=2, col="black") # set at divisions between populations
		  }
		  axis(2, las=2, cex=0.75, pos=0)
		  axis(1, at=c(seq(from = 0.5, to = ncol(plotData + 0.5), by = 1)), labels=sample.names,
			   tick=FALSE, las=2, cex.axis=0.3, adj=0, pos=0, line=0, mgp = c(3, 0.1, 0))
		  axis(1, at=pop.labels, labels=names(table(pop.order.ordered)),
			   tick=FALSE, las=0, cex.axis=0.8, adj=0, pos=0)
		dev.off()
	  } else if(plot.type == "pdf"){
		pdf(file=paste(plot.name, ".apriori.pdf", sep=""), width=plot.width, height=plot.height)
		  par(mgp=c(1.5,1,0), mar = c(4.1, 3.1, 1, 1))
		  barplot(plotData, col=colors, border=border.color,
				  space=0, axes=F, axisname=FALSE, las=2, cex.names=0.75,
				  ylab  = "Assignment Probability", xlab = "Group")
		  for(i in 1:(length(pop.lines)-1)){
			abline(v=pop.lines[i], lwd=2, lty=2, col="black") # set at divisions between populations
		  }
		  axis(2, las=2, cex=0.75, pos=0)
		  axis(1, at=c(seq(from = 0.5, to = ncol(plotData + 0.5), by = 1)), labels=sample.names,
			   tick=FALSE, las=2, cex.axis=0.3, adj=0, pos=0, line=0, mgp = c(3, 0.1, 0))
		  axis(1, at=pop.labels, labels=names(table(pop.order.ordered)),
			   tick=FALSE, las=0, cex.axis=0.8, adj=0, pos=0)
		dev.off()
	  } else {
		png(file=paste(plot.name, ".apriori.png", sep=""), units="in", res=300, width=plot.width, height=plot.height)
		par(mgp=c(1.5,1,0), mar = c(4.1, 3.1, 1, 1))
		barplot(plotData, col=colors, border=border.color,
				space=0, axes=F, axisname=FALSE, las=2, cex.names=0.75,
				ylab  = "Assignment Probability", xlab = "Group")
		for(i in 1:(length(pop.lines)-1)){
		  abline(v=pop.lines[i], lwd=2, lty=2, col="black") # set at divisions between populations
		}
		axis(2, las=2, cex=0.75, pos=0)
		axis(1, at=c(seq(from = 0.5, to = ncol(plotData + 0.5), by = 1)), labels=sample.names,
			 tick=FALSE, las=2, cex.axis=0.3, adj=0, pos=0, line=0, mgp = c(3, 0.1, 0))
		axis(1, at=pop.labels, labels=names(table(pop.order.ordered)),
			 tick=FALSE, las=0, cex.axis=0.8, adj=0, pos=0)
		dev.off()
	  }

	  par(mgp=c(1.5,1,0), mar = c(4.1, 3.1, 1, 1))
	  barplot(plotData, col=colors, border=border.color,
			  space=0, axes=F, axisname=FALSE, las=2, cex.names=0.75,
			  ylab  = "Assignment Probability", xlab = "Group")
	  for(i in 1:(length(pop.lines)-1)){
		abline(v=pop.lines[i], lwd=2, lty=2, col="black") # set at divisions between populations
	  }
	  axis(2, las=2, cex=0.75, pos=0)
	  axis(1, at=c(seq(from = 0.5, to = ncol(plotData + 0.5), by = 1)), labels=sample.names,
		   tick=FALSE, las=2, cex.axis=0.3, adj=0, pos=0, line=0, mgp = c(3, 0.1, 0))
	  axis(1, at=pop.labels, labels=names(table(pop.order.ordered)),
		   tick=FALSE, las=0, cex.axis=0.8, adj=0, pos=0)
	  title()  
  }
}


# function to prep and order DAPC data for plotting, sorts by pop.order
readData <- function(dapc.data, sample.names, pop.order)
{
  data <- cbind(sample.names, pop.order, dapc.data)
  data <- data[order(data[,2]),]

  # calculate number of clusters and taxa
  nclust <- ncol(dapc.data)
  ntax <- nrow(data)

  # build new matrix with only relevant data and labels
  struc_k <- matrix(nrow=nclust, ncol=ntax)
  colnames(struc_k) <- data[,1]
  rownames(struc_k) <- 1:nclust
  for (i in 1:nclust)
  {
    struc_k[i,] <- data[,(2+i)]
  }

  # re-order assignment probs so all plots have similar color schemes
  struc_k_final <- matrix(nrow=nclust, ncol=ntax)
  colnames(struc_k_final) <- data$sample.names
  rownames(struc_k_final) <- 1:nclust
  tracker = NULL
  t = 0
  c = 1
  for(i in 1:ncol(struc_k)){
    m <- max(struc_k[,i])
    for(j in 1:nrow(struc_k)){
      if(m == struc_k[j, i]){
        t = j
      }
    }
    if(!(t %in% tracker)){
      struc_k_final[c,] <- struc_k[t,]
      tracker = c(tracker, t)
      c = c + 1
      #print(tracker)
    }
  }
  for(i in 1:nclust){
    if(!(i %in% tracker)){
      #print(i)
      struc_k_final[c,] <- struc_k[i,]
      c = c + 1
    }
  }
  return(struc_k_final)
}

# function to prep and order DAPC data for plotting, sorts by pop.order then sample.order
readData_sampleOrder <- function(dapc.data, sample.names, pop.order, sample.order)
{
  data <- cbind(sample.names, pop.order, sample.order, dapc.data)
  data <- data[order(data[,2], data[,3]),]

  # calculate number of clusters and taxa
  nclust <- ncol(dapc.data)
  ntax <- nrow(data)

  # build new matrix with only relevant data and labels
  struc_k <- matrix(nrow=nclust, ncol=ntax)
  colnames(struc_k) <- data[,1]
  rownames(struc_k) <- 1:nclust
  for (i in 1:nclust)
  {
    struc_k[i,] <- data[,(3+i)]
  }

  # re-order assignment probs so all plots have similar color schemes
  struc_k_final <- matrix(nrow=nclust, ncol=ntax)
  colnames(struc_k_final) <- data$sample.names
  rownames(struc_k_final) <- 1:nclust
  tracker = NULL
  t = 0
  c = 1
  for(i in 1:ncol(struc_k)){
    m <- max(struc_k[,i])
    for(j in 1:nrow(struc_k)){
      if(m == struc_k[j, i]){
        t = j
      }
    }
    if(!(t %in% tracker)){
      struc_k_final[c,] <- struc_k[t,]
      tracker = c(tracker, t)
      c = c + 1
      #print(tracker)
    }
  }
  for(i in 1:nclust){
    if(!(i %in% tracker)){
      #print(i)
      struc_k_final[c,] <- struc_k[i,]
      c = c + 1
    }
  }
  return(struc_k_final)
}
