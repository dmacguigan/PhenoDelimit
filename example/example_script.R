# example script to run PhenoDelimit

# installation
# install.packages("devtools") # if you need the devtools package
library(devtools)
install_github("dmacguigan/PhenoDelimit", force=TRUE)

library(PhenoDelimit)
library(RColorBrewer)

# step 1: K-means clustering, discriminant analysis, and prep files for CLUMPP
# requires that CLUMPP executable is in the global path
# data were generated using simulate_data.R
wd = "H:/PhenoDelimit/example" # results will be written to new subdirectory "CLUMPP"
data = read.table("H:/PhenoDelimit/example/sim_data.txt", header=TRUE)
n.groups = c(4,2,3,5,5,4)
model.numbers = c(1:6)
models = read.table("H:/PhenoDelimit/example/sim_models.txt", header=TRUE)
perc.var = c(70,80,90)
scale = TRUE
center = TRUE
dapc_clumpp(wd=wd,
            data=data,
            n.groups=n.groups,
            model.numbers=model.numbers,
            models=models,
            perc.var=perc.var,
            scale=scale,
            center=center,
            apriori=FALSE)

# let's run this step again but using a priori population assignments
# we'll come back to these results later
dapc_clumpp(wd=wd,
            data=data,
            n.groups=n.groups,
            model.numbers=model.numbers,
            models=models,
            perc.var=perc.var,
            scale=scale,
            center=center,
            apriori=TRUE)

# step 2: summarize CLUMPP results
# only using our k-means clustering results
wd = "H:/PhenoDelimit/example"
model.numbers = c(1:6)
perc.var = c(70,80,90)

clumpp_results <- read_clumpp_results(wd=wd,
                                      perc.var=perc.var,
                                      model.numbers=model.numbers)

# step 3: plot H' values to compare delimitation models
wd = "H:/PhenoDelimit/example/"
clumpp.data = clumpp_results
colors = c("gray70", "gray30", "black")
plot.type = "png"
plot.name = "H_plot_example"
plot.width = 8
plot.height = 4

plot_clumpp_results(wd=wd,
                    clumpp.data=clumpp.data,
                    colors=colors,
                    plot.name = plot.name,
                    plot.type=plot.type,
                    plot.width=plot.width,
                    plot.height=plot.height)

# step 4: permutation test of significance for H' values
# this may take a while depending on how many models you have
# and how many permutations you wish to perform
wd = "H:/PhenoDelimit/example" # results will be written to new subdirectory "CLUMPP_permuted"
data = read.table("H:/PhenoDelimit/example/sim_data.txt", header=TRUE)
n.groups = c(4,2,3,5,5,4)
model.numbers = c(1:6)
models = read.table("H:/PhenoDelimit/example/sim_models.txt", header=TRUE)
perc.var = c(70,80,90)
scale = TRUE
center = TRUE

dapc_clumpp_permuted(wd=wd,
                     data=data,
                     n.groups=n.groups,
                     model.numbers=model.numbers,
                     models=models,
                     perc.var=perc.var,
                     permutations=100,
                     scale=scale,
                     center=center)

clumpp_perm_df <- read_clumpp_results_permuted(wd=wd,
                                               perc.var=perc.var,
                                               model.numbers=model.numbers,
                                               permutations=100)

H_permutation_plot(wd=wd,
                   clumpp.data=clumpp_results,
                   clumpp.data.permuted=clumpp_perm_df,
                   model.numbers=model.numbers,
                   best.perc.var=70,
                   plot.type="png",
                   plot.prefix="example",
                   plot.width=6,
                   plot.height=4,
                   sig.threshold=0.05)

# step 5: bar plots of discriminant analysis assignment probabilities
models = read.table("H:/PhenoDelimit/example/sim_models.txt", header=TRUE)
wd = "H:/PhenoDelimit/example/"
clumpp.wd = "H:/PhenoDelimit/example/CLUMPP"
sample.names = (1:nrow(models))
sample.plot.groups = models$m1
sample.plot.groups.order = c(1,2,3,4)
sample.order = (nrow(models):1)
best.perc.var = 90
best.model.number = 1
plot.type = "png"
plot.width = 20
plot.height = 4
colors = brewer.pal(n = 5, name = "Set1")
border.color = "gray"

assign_probs_barplot(wd = wd,
                     clumpp.wd = clumpp.wd,
                     sample.names = sample.names,
                     sample.plot.groups = sample.plot.groups,
                     sample.plot.groups.order = sample.plot.groups.order,
                     #sample.order = sample.order,
                     best.perc.var = best.perc.var,
                     best.model.number = best.model.number,
                     plot.type = plot.type,
                     plot.width = plot.width,
                     plot.height = plot.height,
                     colors = colors,
                     border.color = border.color,
                     apriori=FALSE)

# we can also creat a barplot based on results using a priori assignment of individuals to clusters
assign_probs_barplot(wd = wd,
                     clumpp.wd = clumpp.wd,
                     sample.names = sample.names,
                     sample.plot.groups = sample.plot.groups,
                     sample.plot.groups.order = sample.plot.groups.order,
                     #sample.order = sample.order,
                     best.perc.var = best.perc.var,
                     best.model.number = best.model.number,
                     plot.type = plot.type,
                     plot.width = plot.width,
                     plot.height = plot.height,
                     colors = colors,
                     border.color = border.color,
                     apriori=TRUE)

# lets make barplots for all of the other models using a loop
for(m in 2:6){
  best.model.number = m
  sample.names = (1:nrow(models))
  sample.plot.groups = models[,best.model.number]
  sample.plot.groups.order = names(table(sample.plot.groups))
  sample.order = (nrow(models):1)
  best.perc.var = 90

  assign_probs_barplot(wd = wd,
                       clumpp.wd = clumpp.wd,
                       sample.names = sample.names,
                       sample.plot.groups = sample.plot.groups,
                       sample.plot.groups.order = sample.plot.groups.order,
                       sample.order = sample.order,
                       best.perc.var = best.perc.var,
                       best.model.number = best.model.number,
                       plot.type = plot.type,
                       plot.width = plot.width,
                       plot.height = plot.height,
                       colors = colors,
                       border.color = border.color,
                       apriori=FALSE)

  assign_probs_barplot(wd = wd,
                       clumpp.wd = clumpp.wd,
                       sample.names = sample.names,
                       sample.plot.groups = sample.plot.groups,
                       sample.plot.groups.order = sample.plot.groups.order,
                       sample.order = sample.order,
                       best.perc.var = best.perc.var,
                       best.model.number = m,
                       plot.type = plot.type,
                       plot.width = plot.width,
                       plot.height = plot.height,
                       colors = colors,
                       border.color = border.color,
                       apriori=TRUE)
}


# step 6: scatter plot or density plot of discriminant axes
models = read.table("H:/PhenoDelimit/example/sim_models.txt", header=TRUE)

wd = "H:/PhenoDelimit/example/"
clumpp.wd = "H:/PhenoDelimit/example/CLUMPP"
sample.plot.groups = models$m1
sample.plot.groups.order = c(1,2,3,4)
sample.order = (nrow(models):1)
best.perc.var = 90
best.model.number = 1
plot.type = "png"
plot.width = 6
plot.height = 6
colors = brewer.pal(n = 5, name = "Set1")
border.color = "gray"
shapes=c(1:4)

x.axis=1
y.axis=2
plot_discriminant_axes(wd = wd,
                       clumpp.wd = clumpp.wd,
                       sample.plot.groups = sample.plot.groups,
                       sample.plot.groups.order = sample.plot.groups.order,
                       best.perc.var = best.perc.var,
                       best.model.number = best.model.number,
                       plot.type = plot.type,
                       plot.width = plot.width,
                       plot.height = plot.height,
                       colors = colors,
                       shapes = shapes,
                       x.axis = x.axis,
                       y.axis = y.axis,
                       apriori=FALSE)

y.axis=3
plot_discriminant_axes(wd = wd,
                       clumpp.wd = clumpp.wd,
                       sample.plot.groups = sample.plot.groups,
                       sample.plot.groups.order = sample.plot.groups.order,
                       best.perc.var = best.perc.var,
                       best.model.number = best.model.number,
                       plot.type = plot.type,
                       plot.width = plot.width,
                       plot.height = plot.height,
                       colors = colors,
                       shapes = shapes,
                       x.axis = x.axis,
                       y.axis = y.axis,
                       apriori=FALSE)

# now make the same plots, but using the a priori group assignments
x.axis=1
y.axis=2
plot_discriminant_axes(wd = wd,
                       clumpp.wd = clumpp.wd,
                       sample.plot.groups = sample.plot.groups,
                       sample.plot.groups.order = sample.plot.groups.order,
                       best.perc.var = best.perc.var,
                       best.model.number = best.model.number,
                       plot.type = plot.type,
                       plot.width = plot.width,
                       plot.height = plot.height,
                       colors = colors,
                       shapes = shapes,
                       x.axis = x.axis,
                       y.axis = y.axis,
                       apriori=TRUE)

y.axis=3
plot_discriminant_axes(wd = wd,
                       clumpp.wd = clumpp.wd,
                       sample.plot.groups = sample.plot.groups,
                       sample.plot.groups.order = sample.plot.groups.order,
                       best.perc.var = best.perc.var,
                       best.model.number = best.model.number,
                       plot.type = plot.type,
                       plot.width = plot.width,
                       plot.height = plot.height,
                       colors = colors,
                       shapes = shapes,
                       x.axis = x.axis,
                       y.axis = y.axis,
                       apriori=TRUE)

# step 7: plot discriminant axis loadings and write to table
wd = "H:/PhenoDelimit/example/"
clumpp.wd = "H:/PhenoDelimit/example/CLUMPP"
best.perc.var = 90
best.model.number = 1
plot.type = "png"
plot.width = 6
plot.height = 6

axis=1
discriminant_loading(wd = wd,
                     clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var,
                     best.model.number = best.model.number,
                     plot.type = plot.type,
                     plot.width = plot.width,
                     plot.height = plot.height,
                     axis = axis,
                     apriori=FALSE)

axis=2
discriminant_loading(wd = wd,
                     clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var,
                     best.model.number = best.model.number,
                     plot.type = plot.type,
                     plot.width = plot.width,
                     plot.height = plot.height,
                     axis = axis,
                     apriori=FALSE)

axis=3
discriminant_loading(wd = wd,
                     clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var,
                     best.model.number = best.model.number,
                     plot.type = plot.type,
                     plot.width = plot.width,
                     plot.height = plot.height,
                     axis = axis,
                     apriori=FALSE)


# now make the same plots, but using the a priori group assignments
axis=1
discriminant_loading(wd = wd,
                     clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var,
                     best.model.number = best.model.number,
                     plot.type = plot.type,
                     plot.width = plot.width,
                     plot.height = plot.height,
                     axis = axis,
                     apriori=TRUE)

axis=2
discriminant_loading(wd = wd,
                     clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var,
                     best.model.number = best.model.number,
                     plot.type = plot.type,
                     plot.width = plot.width,
                     plot.height = plot.height,
                     axis = axis,
                     apriori=TRUE)

axis=3
discriminant_loading(wd = wd,
                     clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var,
                     best.model.number = best.model.number,
                     plot.type = plot.type,
                     plot.width = plot.width,
                     plot.height = plot.height,
                     axis = axis,
                     apriori=TRUE)


