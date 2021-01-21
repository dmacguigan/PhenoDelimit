# discriminatnt analysis of principal components for Elep meristic data
library(adegenet)
library(ggplot2)
library(reshape2)
library(svglite)
library(plyr)
library(RColorBrewer)
library(ggpubr)

wd = "H:/DJM_NearLabMac/dmacguigan/Documents/NearLab/LepidumProject/DAPC_meristic_CLUMPP"
setwd(wd)

data <- read.csv("Elep_meristicsFINAL.csv", header=TRUE)
data <- data[,-33] # remove SL variable
data <- data[,-28:-32] # remove squamation variables
data <- data[,-24] # remove anal spines variable
data <- data[,-18] # remove second dorsal-anal fin variable
data$caudalPeduncleScales <- data$CaudalPeduncleAboveL.L. + data$CaudalPeduncleBelowL.L. - 1 # combine scales above and below caudal peduncle
data <- data[,-18:-19] # remove scales above and below LL variables

# drop Medina (aka San Antonio) samples because of low sample size
data <- data[-which(data$Clade == "G_Medina"),]

# get data only
d <- data[,16:ncol(data)]

# remove rows with NAs
d.complete <- na.omit(d)
data.complete <- data[!(rowSums(is.na(d)) > 0),]

# perform a box-cox transformation of environmental data
#trans <- preProcess(d.complete, method=c("BoxCox", "center", "scale"))
#PC <- predict(trans, d.complete)

# quick PCA of environmental variables
#meristic.PCA.boxCox <- prcomp(x=PC)
#meristic.PCA <- prcomp(x=d.complete, scale=TRUE, center=TRUE)

#autoplot(meristic.PCA.boxCox, data = data.complete, colour = 'Species', loadings=TRUE, loadings.label=TRUE)
#autoplot(meristic.PCA, data = data.complete, colour = 'Species', loadings=TRUE, loadings.label=TRUE)


########################################################################
# function to create CLUMPP indfiles for DAPC using model assignment vs K-means clusters DAPC probs
clumpp_prep <- function(data, modelNumber, k, model, perc.pca){
  # make sure model variable is OK
  model <- as.numeric(as.factor(model))
  # model DAPC
  dapc.model <- dapc(data, model, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.pca, n.da=10000, center=TRUE, scale=TRUE)
  # K-means DAPC
  kmeans.cluster <- find.clusters(data, max.n.clust = 30, n.start = 1000, n.iter=1e6, n.pca=10000, n.clust=k, center=TRUE, scale=TRUE) #retain all PCs
  dapc.kmeans <- dapc(data, kmeans.cluster$grp, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.pca, n.da=10000, center=TRUE, scale=TRUE)

  # create CLUMPP indfile
  d1 <- data.frame(1:nrow(dapc.model$posterior),
                   1:nrow(dapc.model$posterior),
                   rep("(x)", nrow(dapc.model$posterior)),
                   model,
                   rep(":", nrow(dapc.model$posterior)),
                   dapc.model$posterior)
  d2 <- data.frame(1:nrow(dapc.kmeans$posterior),
                   1:nrow(dapc.kmeans$posterior),
                   rep("(x)", nrow(dapc.kmeans$posterior)),
                   model,
                   rep(":", nrow(dapc.kmeans$posterior)),
                   dapc.kmeans$posterior)
  colnames(d2) <- colnames(d1)
  d <- rbind(d1, d2)
  write.table(file = paste("m", modelNumber,"_percPCA-", perc.pca, ".indfile", sep=""), x=d,
              quote=FALSE, col.names = FALSE, row.names=FALSE)
}

# function to create CLUMPP indfiles for model assignment vs model DAPC probs
clumpp_prep_model_vs_pop <- function(data, modelNumber, k, model, perc.pca){
  # make sure model variable is OK
  model <- as.numeric(as.factor(model))
  # model DAPC
  dapc.model <- dapc(data, model, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.pca, n.da=10000, center=TRUE, scale=TRUE)
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
  d2 <- data.frame(1:nrow(dapc.model$posterior),
                   1:nrow(dapc.model$posterior),
                   rep("(x)", nrow(dapc.model$posterior)),
                   model,
                   rep(":", nrow(dapc.model$posterior)),
                   dapc.model$posterior)
  colnames(d2) <- colnames(d1)
  d <- rbind(d1, d2)
  write.table(file = paste("m", modelNumber, "_percPCA-", perc.pca, "_modelVsPop.indfile", sep=""), x=d,
              quote=FALSE, col.names = FALSE, row.names=FALSE)
}

# function to create CLUMPP indfiles for model assignment vs K-means clusters DAPC probs
clumpp_prep_kmeans_vs_pop <- function(data, modelNumber, k, model, perc.pca){
  # make sure model variable is OK
  model <- as.numeric(as.factor(model))
  # K-means DAPC
  kmeans.cluster <- find.clusters(data, max.n.clust = 30, n.start = 1000, n.iter=1e6, n.pca=10000, n.clust=k, center=TRUE, scale=TRUE) #retain all PCs
  dapc.kmeans <- dapc(data, kmeans.cluster$grp, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.pca, n.da=10000, center=TRUE, scale=TRUE)
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
  write.table(file = paste("m", modelNumber, "_percPCA-", perc.pca, "_kmeansVsPop.indfile", sep=""), x=d,
              quote=FALSE, col.names = FALSE, row.names=FALSE)
}


########################################################################
# create CLUMPP indfiles
setwd(wd)
setwd("./CLUMPP")
# model info
nSpecies <- c(2,4,8,3,4,7) # number of species in each model, slightly different due to dropping San Antionio
modelNumber <- c(2,3,4,5,6,7) # model numbers
models <- data.complete[,10:15] # delimitation models
perc.pca <- c(70,80,90) # percentage of varaince to retain
for(i in 1:length(nSpecies)){
  for(j in 1:length(perc.pca)){
    clumpp_prep(data=d.complete, modelNumber=modelNumber[i], k=nSpecies[i], model=as.factor(models[,i]), perc.pca=perc.pca[j])
    clumpp_prep_model_vs_pop(data=d.complete, modelNumber=modelNumber[i], k=nSpecies[i], model=as.factor(models[,i]), perc.pca=perc.pca[j])
    clumpp_prep_kmeans_vs_pop(data=d.complete, modelNumber=modelNumber[i], k=nSpecies[i], model=as.factor(models[,i]), perc.pca=perc.pca[j])
  }
}

########################################################################
# read in results from miscfiles
setwd(wd)
setwd("./CLUMPP")
suffix=c("kmeansVsPop", "modelVsPop")
c1=numeric()
c2=numeric()
c3=character()
c4=numeric()
for(x in 1:length(suffix)){
  for(i in 1:length(perc.pca)){
    for(j in 1:length(modelNumber)){
      fileName <- Sys.glob(paste("m", modelNumber[j], "_percPCA-", perc.pca[i], "_", suffix[x], ".miscfile", sep=""))
      temp <- grep("The highest value of H' is:", readLines(fileName), value=TRUE)
      H <- as.numeric(strsplit(temp, " ")[[1]][7]) # get only H' number
      c1 <- c(c1, modelNumber[j])
      c2 <- c(c2, perc.pca[i])
      c3 <- c(c3, suffix[x])
      c4 <- c(c4, H)
    }
  }
}
H_df <- data.frame(model=c1, perc.pca=c2, comparison=c3, plotGroup=paste(c2, c3, sep="_"), H=c4)

# which model is the max for each plotting group?
maxes <- aggregate(H_df$H, by = list(H_df$plotGroup), max)
H_df$shape <- ifelse((H_df$H %in% maxes$x), 16, 1)

########################################################################
# plot all results
p <- ggplot(data = H_df, aes(y=H, x=model, pch=shape, col=comparison)) +
  scale_shape_identity() +
  geom_point(size=5) +
  geom_line(aes(linetype=as.factor(perc.pca)), lwd=1.2) +
  scale_color_manual(values = c("gray", "black"), labels=c("K-means", "model"), name="cluster\nassignment") +
  scale_linetype_manual(name="percent\nvariance\nretained", values=c(1, 2, 3)) +
  ylab(label = "H'") +
  ylim(c(0,1)) +
  theme_minimal()
  #theme(legend.position = "none")
svglite(file = "Elep_meristic_HPrime_all.svg", width = 7.5, height = 4)
print(p)
dev.off()


# plot only max percent variance retained
H_df_90 <- subset(H_df, perc.pca == 90)
p <- ggplot(data = H_df_90, aes(y=H, x=model, pch=shape, col=comparison)) +
  scale_shape_identity() +
  geom_point(size=5) +
  geom_line() +
  scale_color_manual(values = c("gray", "black"), labels=c("K-means", "model"), name="cluster\nassignment") +
  ylab(label = "H'") +
  ylim(c(0,1)) +
  theme_minimal()
#theme(legend.position = "none")
svglite(file = "Elep_meristic_HPrime_90PercVar.svg", width = 7.5, height = 1.5)
print(p)
dev.off()


########################################################################
# make heatmap plots of best result
# K-means clustering

# function to make heatmap plot of assignment probabilties CLUMPP
# use K-means for DAPC a prior cluster assignment
plot_DAPC <- function(data, y.axis, k, model, perc.pca){
  # make sure model variable is OK
  model <- as.numeric(as.factor(model))
  # K-means DAPC
  kmeans.cluster <- find.clusters(data, max.n.clust = 30, n.start = 1000, n.iter=1e6, n.pca=10000, n.clust=k, center=TRUE, scale=TRUE) #retain all PCs
  dapc.kmeans <- dapc(data, kmeans.cluster$grp, var.conrib=TRUE, var.loadings=TRUE, perc.pca=perc.pca, n.da=10000, center=TRUE, scale=TRUE)

  h <- data.frame(dapc.kmeans$posterior)
  h$Clade <- y.axis
  h <- aggregate(h[,1:k], list(h$Clade), mean)
  h <- reshape2::melt(h)
  colnames(h) <- c("River", "Species", "probability")

  hm <- ggplot(h, aes(x=Species, y=River)) +
    geom_tile(aes(fill=probability)) +
    geom_text(aes(label = round(probability, 2))) +
    scale_fill_gradient(low="#800026", high="#ffffcc", breaks=seq(0,1, 0.25), limits=c(0,1)) +
    theme_minimal()

  svglite(file="Elep_meristic_hm.svg", width=5, height=5)
  print(hm)
  dev.off()

  return(dapc.kmeans)

}

models <- data.complete[,10:15] # delimitation models
data=d.complete
k=3
model=as.factor(models[,4])
perc.pca=90
y.axis=factor(data.complete$Clade)
y.axis <- revalue(y.axis, c("A_Pecos"="Pecos", "B_Concho"="Concho",
                          "C_SanSaba"="San Saba", "D_Llano"="Llano",
                          "E_Pedernales"="Pedernales", "F_Guadalupe"="Guadalupe",
                          "H_Frio"="Frio", "I_Nueces"="Nueces"))
y.axis <- factor(y.axis, levels=c("Frio", "Nueces", "Guadalupe", "Pedernales",
                               "Llano", "San Saba", "Concho", "Pecos"))

setwd(wd)
dapc.kmeans <- plot_DAPC(data=data, y.axis=y.axis, k=k, model=model, perc.pca=perc.pca)


########################################################################
# make barplot of best result posterior probabilities

# function to prep and order DAPC data for plotting
readData <- function(dapc.data, sample.data, order.col, name.col)
{
  data <- cbind(sample.data,dapc.data)
  data <- data[order(data[,order.col]),]

  data

  # calculate number of clusters and taxa
  nclust <- (ncol(data)-ncol(sample.data))
  ntax <- nrow(data)

  # build new matrix with only relevant data and labels
  struc_k <- matrix(nrow=nclust, ncol=ntax)
  colnames(struc_k) <- data[,1]
  rownames(struc_k) <- 1:nclust
  for (i in 1:nclust)
  {
    struc_k[i,] <- data[,(ncol(sample.data)+i)]
  }

  # re-order assignment probs so all plots have similar color schemes
  struc_k_final <- matrix(nrow=nclust, ncol=ntax)
  colnames(struc_k_final) <- data[,name.col]
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
      print(tracker)
    }
  }
  for(i in 1:nclust){
    if(!(i %in% tracker)){
      print(i)
      struc_k_final[c,] <- struc_k[i,]
      c = c + 1
    }
  }
  struc_k_final
}

# set up order of rivers for barplot
order.col=7
name.col=7
sample.data=data.complete[,1:8]
sample.data[,order.col]=factor(sample.data[,order.col])
sample.data[,order.col] <- revalue(sample.data[,order.col], c("A_Pecos"="Pecos", "B_Concho"="Concho",
                            "C_SanSaba"="San Saba", "D_Llano"="Llano",
                            "E_Pedernales"="Pedernales", "F_Guadalupe"="Guadalupe",
                            "H_Frio"="Frio", "I_Nueces"="Nueces"))
sample.data[,order.col] <- factor(sample.data[,order.col], levels=rev(c("Frio", "Nueces", "Guadalupe", "Pedernales",
                                  "Llano", "San Saba", "Concho", "Pecos")))

# prep data for plotting
barplot.data <- readData(dapc.data=dapc.kmeans$posterior, sample.data=sample.data, order.col=order.col, name.col=name.col)

# get positions of breaks between populations and species
pop.col=7
sp.col=8
sample.data.ordered <- sample.data[order(sample.data[,pop.col]),]
sample.data.ordered[,sp.col]=factor(sample.data.ordered[,sp.col])
sample.data.ordered[,sp.col] <- revalue(sample.data.ordered[,sp.col],
                                            c("A_Pecos"="Pecos", "B_Concho-SanSaba"="Upper Colorado", "C_EastEdPlat"="Eastern Edwards Plateau"))
sample.data.ordered[,sp.col] <- factor(sample.data.ordered[,sp.col], levels=c("Pecos", "Upper Colorado", "Eastern Edwards Plateau"))
pop.lines <- cumsum(table(sample.data.ordered[,pop.col]))
pop.labels <- (pop.lines - table(sample.data.ordered[,pop.col])) + (table(sample.data.ordered[,pop.col])/2)
sp.lines <- cumsum(table(sample.data.ordered[,sp.col]))
sp.labels <- (sp.lines - table(sample.data.ordered[,sp.col])) + (table(sample.data.ordered[,sp.col])/2)

setwd(wd)
svglite(file="Elep_meristic_barplot.svg", width=8, height=3)
#pdf(file="Elep_meristic_barplot.pdf", width=10, height=3)
barplot(barplot.data, col=c("#fc8d62", "#66c2a5", "#8da0cb"), border=NA,
        space=0, axes=F, axisname=FALSE, las=2, cex.names=0.75)
for(i in 1:(length(pop.lines)-1)){
  abline(v=pop.lines[i], lwd=1, lty=2, col="black") # set at divisions between populations
}
for(i in 1:(length(sp.lines)-1)){
  abline(v=sp.lines[i], lwd=2, lty=1, col="black") # set at divisions between species
}
axis(2, las=2, cex=0.75, pos=0)
axis(1, at=pop.labels,
     labels=c("", "Concho", "San Saba", "Llano", "Pedernales", "Guadalupe", "Frio", "Nueces"),
     tick=FALSE, las=0, cex.axis=0.5, adj=0, pos=0, padj=-1.5)
axis(1, at=sp.labels,
     labels=c("Pecos", "Upper Colorado", "Eastern Edwards Plateau"),
     tick=FALSE, las=0, cex.axis=1, adj=0, pos=0, padj=1)
dev.off()

########################################################################
# create scatter plot of LD axes

# function to return LD 1 vs 2 biplot for 3 species model
LD.biplot <- function(data, speciesColors, species){
  clusterShapes <- c(0,1,2)
  p <- data.frame(data$ind.coord)
  p$species <- species

  find_hull <- function(df) df[chull(df$LD1, df$LD2), ]
  hulls <- ddply(p, "species", find_hull)

  p$cluster <- data$grp

  plot <- ggplot() +
    geom_polygon(data=hulls, aes(x=LD1, y=LD2, color=species, fill=species), alpha=0.2) +
    geom_point(data=p, aes(x=LD1, y=LD2, color=species, fill=species, shape=cluster), cex=3, alpha=0.8) +
    scale_color_manual(values = speciesColors, labels = c("Pecos", "Upper Colorado", "Eastern Edwards\nPlateau"), name = "Species") +
    scale_fill_manual(values = speciesColors, guide="none") +
    scale_shape_manual(values = clusterShapes, labels = c("1","2","3"), name = "K-means\nCluster") +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), legend.box = "horizontal") +
    theme_minimal() +
    xlab("Discriminant axis 1") +
    ylab("Discriminant axis 2")

  plot

  return(plot)
}

setwd(wd)

sp <- revalue(data.complete$Species,  c("A_Pecos"="Pecos", "B_Concho-SanSaba"="Upper Colorado", "C_EastEdPlat"="Eastern Edwards Plateau"))
sp <- factor(sp, levels=c("Pecos", "Upper Colorado", "Eastern Edwards Plateau"))

p <- LD.biplot(data = dapc.kmeans,
               speciesColors = c("#66c2a5", "#8da0cb", "#fc8d62"),
               species = sp)
p
svglite(file = "Elep_meristic_scatter.svg", width = 7.5, height = 4)
print(p)
dev.off()

# function to plot discriminant axis loadings
# and write to files
plot.DA.loading <- function(data, axis, name){
  pdf(file=paste(name, ".loading.DA", axis, ".pdf", sep=""), width=8, height = 6)
  loadingplot(data$var.contr, axis=axis)
  dev.off()
  d <- data$var.contr[,axis]
  d2 <- data$var.load[,axis]
  t <- cbind(d,d2)
  t <- as.data.frame(t)
  t <- t[order(t$d, decreasing = TRUE),]
  best <- t[1:10,]
  colnames(best) <- c("contribution", "loading")
  write.table(x = best, file = paste(name, ".loading.DA", axis, ".txt", sep=""), quote=FALSE)
}

setwd(wd)
plot.DA.loading(dapc.kmeans, 1, name="dapc.90")
plot.DA.loading(dapc.kmeans, 2, name="dapc.90")


