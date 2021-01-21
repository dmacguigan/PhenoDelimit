# script to simulate data
# simulates K groups with roughly equal sample sizes

library(ggplot2)

# working directory
wd = "H:/NearLab/PhenoDelimit/example"

# how many clusters to generate?
true_clust = 4
# how many samples total?
samples = 300
# how many variables to generate?
nVar = 20
# normal distribution sd range
var_sd_min = 1
var_sd_max = 10
# variable normal distribution mean min and max
var_mean_min = 5
var_mean_max = 30

# generate variable normal distribution means
temp <- list()
for(i in 1:true_clust){
  temp[[i]] <- runif(nVar, var_mean_min, var_mean_max) # choose random value from uniform distribution between var_mean_min and var_mean_max
}
var_means <- as.data.frame(do.call(rbind, temp))

# generate variable normal distribution standard deviations, use same sd for each variable
var_sds <- runif(nVar, var_sd_min, var_sd_max) # choose random value from uniform distribution between var_sd_min and var_sd_max

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

# quick plot of one variable
plot(sim_data$cluster, sim_data$var2)

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

