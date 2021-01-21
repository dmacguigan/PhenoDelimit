# PhenoDelimit
R package to create "structure-like" barplots and compare species delimitation models with multivariate phenotype data

This package requires an external program, CLUMPP, which cana be dowloaded here: https://rosenberglab.stanford.edu/clumpp.html


To download and install this package
```
install.packages("devtools")
library(devtools)
install_github("dmacguigan/PhenoDelimit")
library(PhenoDelimit)
```

Example script for how to run the package. Example files in "example" directory, so clone the repository first.
```
# example script to run PhenoDelimit

# step 1: K-means clustering, discriminant analysis, and prep files for CLUMPP
wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
data = read.table("H:/NearLab/PhenoDelimit/example/sim_data.txt", header=TRUE)
n.groups = c(4,2,3,5,5,4)
model.numbers = c(1:6)
models = read.table("H:/NearLab/PhenoDelimit/example/sim_models.txt", header=TRUE)
perc.var = c(70,80,90)
scale = TRUE
center = TRUE

clumpp_prep(wd=wd, data=data, n.groups=n.groups, model.numbers=model.numbers, models=models, perc.var=perc.var, scale=scale, center=center)

# step 2: run CLUMPP, must be done externally
# can use run_CLUMPP.sh in "PhenoDelimit/example/CLUMPP" to run CLUMPP on all files in a directory

# step 3: summarize CLUMPP
wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
model.numbers = c(1:6)
perc.var = c(70,80,90)

clumpp_results <- read_clumpp_results(wd=wd, perc.var=perc.var, model.numbers=model.numbers)

# step 4: plot H' values to compare delimitation models
wd = "H:/NearLab/PhenoDelimit/example/"
clumpp.data = clumpp_results
colors = brewer.pal(n = 3, name = "Set1")
plot.type = "pdf"
plot.name = "H_plot_example"
plot.width = 8
plot.height = 4

plot_clumpp_results(wd=wd, clumpp.data=clumpp.data, colors=colors, plot.name = plot.name, plot.type=plot.type, plot.width=plot.width, plot.height=plot.height)

# step 5: barplot
models = read.table("H:/NearLab/PhenoDelimit/example/sim_models.txt", header=TRUE)

wd = "H:/NearLab/PhenoDelimit/example/"
clumpp.wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
sample.names = (1:nrow(models))
sample.plot.groups = models$m1
sample.plot.groups.order = c(1,2,3,4)
sample.order = (nrow(models):1)
best.perc.var = 90
best.model.number = 1
plot.type = "pdf"
plot.width = 8
plot.height = 6
plot.name = "barplot_example"
colors = brewer.pal(n = 5, name = "Set1")
border.color = "gray"

assignProbs_barplot(wd = wd, clumpp.wd = clumpp.wd, sample.names = sample.names,
                    sample.plot.groups = sample.plot.groups, sample.plot.groups.order = sample.plot.groups.order,
                    #sample.order = sample.order,
                    best.perc.var = best.perc.var, best.model.number = best.model.number,
                    plot.type = plot.type, plot.width = plot.width, plot.height = plot.height,
                    plot.name = plot.name, colors = colors, border.color = border.color)
```
