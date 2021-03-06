# script to simulate data
# K groups with roughly equal sample sizes
# generate 20 variables
# var 1-10 are allowed to have large group differences, small sds
# var 11-20 have small group differences, large sds

library(ggplot2)
library(reshape2)
library(RColorBrewer)

# working directory
wd = "H:/NearLab/PhenoDelimit/example"

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


# basic PCA of simulated data, show PC1 and PC2
sim_pca <- prcomp(sim_data[,2:(nVar+1)], center=TRUE, scale=TRUE)
screeplot(sim_pca)
plot_data <- cbind(cluster=as.factor(sim_data$cluster), as.data.frame(sim_pca$x))
p <- ggplot(data=plot_data, aes(x=PC1, y=PC2, color=cluster)) +
  geom_point()
p

# write simulated data to file
setwd(wd)
write.table(sim_data[,2:(nVar+1)], file = "sim_data.txt", quote=FALSE, row.names = FALSE)

# create a few "delimitation models" to test
# ONLY WORKS FOR true_clust=4

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
setwd(wd)
write.table(models, file = "sim_models.txt", quote=FALSE, row.names = FALSE)

# read data if already created
setwd(wd)
sim_data <- read.table(file = "sim_data.txt", header=TRUE)
sim_data$grp <- as.factor(read.table(file = "sim_models.txt", header=TRUE)[,1])

sim_data_melt <- melt(sim_data, id.vars = "grp")

# boxplot of all variables
p <- ggplot(sim_data_melt, aes(x=grp, y=value, group=grp)) +
  geom_boxplot(aes(fill=grp)) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set1"), name="Group") +
  facet_wrap(~ variable, nrow=2, ncol=10) +
  theme_light()

png(filename = "sim_data_boxplots.png", width=10, height=6, units = "in", res = 300)
p
dev.off()


