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

Example script. Data are in "example" directory, so clone the repository first.

## step 1: K-means clustering, discriminant analysis, and prep files for CLUMPP

```
# data were generated using simulate_data.R in the example folder
library(RColorBrewer)

wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
data = read.table("H:/NearLab/PhenoDelimit/example/sim_data.txt", header=TRUE)
n.groups = c(4,2,3,5,5,4)
model.numbers = c(1:6)
models = read.table("H:/NearLab/PhenoDelimit/example/sim_models.txt", header=TRUE)
perc.var = c(70,80,90)
scale = TRUE
center = TRUE

clumpp_prep(wd=wd, data=data, n.groups=n.groups, model.numbers=model.numbers, models=models, perc.var=perc.var, scale=scale, center=center)
```
## step 2: run CLUMPP, must be done externally
Can use run_CLUMPP.sh in "PhenoDelimit/example/CLUMPP" to run CLUMPP on all files in a directory

## step 3: summarize CLUMPP
```
wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
model.numbers = c(1:6)
perc.var = c(70,80,90)

clumpp_results <- read_clumpp_results(wd=wd, perc.var=perc.var, model.numbers=model.numbers)
```
## step 4: plot H' values to compare delimitation models
```
wd = "H:/NearLab/PhenoDelimit/example/"
clumpp.data = clumpp_results
colors = brewer.pal(n = 3, name = "Set1")
plot.type = "png"
plot.name = "H_plot_example"
plot.width = 8
plot.height = 4

plot_clumpp_results(wd=wd, clumpp.data=clumpp.data, colors=colors, plot.name = plot.name, plot.type=plot.type, plot.width=plot.width, plot.height=plot.height)
```
![H_plot_example](/example/H_plot_example.png)

## step 5: bar plots of discriminant analysis assignment probabilities
```
models = read.table("H:/NearLab/PhenoDelimit/example/sim_models.txt", header=TRUE)

wd = "H:/NearLab/PhenoDelimit/example/"
clumpp.wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
sample.names = (1:nrow(models))
sample.plot.groups = models$m1
sample.plot.groups.order = c(1,2,3,4)
sample.order = (nrow(models):1)
best.perc.var = 90
best.model.number = 1
plot.type = "png"
plot.width = 20
plot.height = 6
plot.name = "barplot_example"
colors = brewer.pal(n = 5, name = "Set1")
border.color = "gray"

assign_probs_barplot(wd = wd, clumpp.wd = clumpp.wd, sample.names = sample.names,
                    sample.plot.groups = sample.plot.groups, sample.plot.groups.order = sample.plot.groups.order,
                    #sample.order = sample.order,
                    best.perc.var = best.perc.var, best.model.number = best.model.number,
                    plot.type = plot.type, plot.width = plot.width, plot.height = plot.height,
                    plot.name = plot.name, colors = colors, border.color = border.color)
```
![barplot_example](/example/barplot_example.png)

## step 6: scatter plot or density plot of discriminant axes
```
models = read.table("H:/NearLab/PhenoDelimit/example/sim_models.txt", header=TRUE)

wd = "H:/NearLab/PhenoDelimit/example/"
clumpp.wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
sample.plot.groups = models$m1
sample.plot.groups.order = c(1,2,3,4)
sample.order = (nrow(models):1)
best.perc.var = 90
best.model.number = 1
plot.type = "png"
plot.width = 6
plot.height = 6
plot.name = "scatter_plot_example"
colors = brewer.pal(n = 5, name = "Set1")
border.color = "gray"
shapes=c(1:4)
x.axis=1
y.axis=2

plot_discriminant_axes(wd = wd, clumpp.wd = clumpp.wd,
                       sample.plot.groups = sample.plot.groups, sample.plot.groups.order = sample.plot.groups.order,
                       best.perc.var = best.perc.var, best.model.number = best.model.number,
                       plot.type = plot.type, plot.width = plot.width, plot.height = plot.height,
                       plot.name = plot.name, colors = colors, shapes = shapes,
                       x.axis = x.axis, y.axis = y.axis)

```
![scatter_plot_example](/example/scatter_plot_example.png)


## step 7: plot discriminant axis loadings and write table of variable contributions and loadings
```
wd = "H:/NearLab/PhenoDelimit/example/"
clumpp.wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
best.perc.var = 90
best.model.number = 1
plot.type = "png"
plot.width = 6
plot.height = 6
axis=1

discriminant_loading(wd = wd, clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var, best.model.number = best.model.number,
                     plot.type = plot.type, plot.width = plot.width, plot.height = plot.height, axis = axis)
```
![DA_1_loading](/example/DA_1_loading.png)




