# PhenoDelimit
R package to create "STRUCTURE-esque" barplots and compare species delimitation models with multivariate phenotype data.

This package requires an external program, CLUMPP, which can be [dowloaded here](https://rosenberglab.stanford.edu/clumpp.html)

To download and install this package
```
install.packages("devtools")
library(devtools)
install_github("dmacguigan/PhenoDelimit")
library(PhenoDelimit)
```

## Example analysis
Requires that CLUMPP executable is in the global path.

Data are in "example" directory, so clone the repository first. Data were generated using simulate_data.R in the example folder
Dataset consists 20 continuous variables for 300 idividuals divided into four groups.
Data were simulated by drawing from normal disrtibutions.
For each variable, means were allowed to vary between groups, but standard deviations were kept constant.

Variables 1-10 had a large range of group means (between 5 and 30) and small standard deviations (between 0.5 and 1). These variables are expected to be more informative about group membership.

Variables 11-20 had a small range of group means (between 18 and 22) and large standard deviations (between 3 and 10). Thus, these variables are expected to be less informative about group membership.
<br/>
<br/>
### step 1: K-means clustering, discriminant analysis of principal components, and run CLUMPP
Use CLUMPP to compare discriminant analysis assignment probabilities to each delimitation model. CLUMPP calculates `H'`, a matrix similarity metric bounded between 0 and 1.
A value of 1 indicates a perfect match between two matrices.
See [this paper in Bioinformatics](https://academic.oup.com/bioinformatics/article/23/14/1801/188285) for more info about the CLUMPP algorithm.

In this analysis, we provide CLUMPP with two matrices.
First, a matrix where each individual is assigned with probability 1.0 to a group based on the delimitation model.
Second, a matrix of individual assignment probabilities fron discriminant analysis using groups defined by K-means clustering.
By comparing the two matrices with CLUMPP, we can determine which delimitation model is the best match to K-means clusters, indicative of real group differences in the data.

```
library(RColorBrewer)

wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
data = read.table("H:/NearLab/PhenoDelimit/example/sim_data.txt", header=TRUE)
n.groups = c(4,2,3,5,5,4)
model.numbers = c(1:6)
models = read.table("H:/NearLab/PhenoDelimit/example/sim_models.txt", header=TRUE)
perc.var = c(70,80,90)
scale = TRUE
center = TRUE

dapc_clumpp(wd=wd, data=data, n.groups=n.groups, model.numbers=model.numbers, models=models, perc.var=perc.var, scale=scale, center=center)
```
<br/>
<br/>
### step 2: summarize CLUMPP
```
wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
model.numbers = c(1:6)
perc.var = c(70,80,90)

clumpp_results <- read_clumpp_results(wd=wd, perc.var=perc.var, model.numbers=model.numbers)
```
<br/>
<br/>
### step 3: plot H' values to compare delimitation models
```
wd = "H:/NearLab/PhenoDelimit/example/"
clumpp.data = clumpp_results # from previous step
colors = brewer.pal(n = 3, name = "Set1")
plot.type = "png"
plot.name = "H_plot_example"
plot.width = 8
plot.height = 4

plot_clumpp_results(wd=wd, clumpp.data=clumpp.data, colors=colors, plot.name = plot.name, plot.type=plot.type, plot.width=plot.width, plot.height=plot.height)
```
![H_plot_example](/example/H_plot_example.png)

For this simulated datset, model 1 is the "true" model which generated the data.
Models 2-5 are tweaked versions of model 1 with groups merged or split.
Model 6 randomly shuffled the "true" group assignments from model 1.
<br/>
<br/>
### step 4: bar plots of discriminant analysis assignment probabilities
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
<br/>
<br/>
### step 5: scatter plot or density plot of discriminant axes
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
colors = brewer.pal(n = 5, name = "Set1")
border.color = "gray"
shapes=c(1:4)


# Disriminant axis 1 vs 2
x.axis=1
y.axis=2

plot_discriminant_axes(wd = wd, clumpp.wd = clumpp.wd,
                       sample.plot.groups = sample.plot.groups, sample.plot.groups.order = sample.plot.groups.order,
                       best.perc.var = best.perc.var, best.model.number = best.model.number,
                       plot.type = plot.type, plot.width = plot.width, plot.height = plot.height,
                       colors = colors, shapes = shapes,
                       x.axis = x.axis, y.axis = y.axis)
```
![DA1-DA2_scatter.png](/example/DA1-DA2_scatter.png)

```
# Disriminant axis 1 vs 3
x.axis=1
y.axis=3

plot_discriminant_axes(wd = wd, clumpp.wd = clumpp.wd,
                       sample.plot.groups = sample.plot.groups, sample.plot.groups.order = sample.plot.groups.order,
                       best.perc.var = best.perc.var, best.model.number = best.model.number,
                       plot.type = plot.type, plot.width = plot.width, plot.height = plot.height,
                       colors = colors, shapes = shapes,
                       x.axis = x.axis, y.axis = y.axis)
```
![DA1-DA3_scatter.png](/example/DA1-DA3_scatter.png)
<br/>
<br/>
### step 6: plot discriminant axis loadings and write table of variable contributions and loadings
```
wd = "H:/NearLab/PhenoDelimit/example/"
clumpp.wd = "H:/NearLab/PhenoDelimit/example/CLUMPP"
best.perc.var = 90
best.model.number = 1
plot.type = "png"
plot.width = 6
plot.height = 6

# Variable loading on discriminant axis 1
axis=1

discriminant_loading(wd = wd, clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var, best.model.number = best.model.number,
                     plot.type = plot.type, plot.width = plot.width, plot.height = plot.height, axis = axis)
```
![DA_1_loading](/example/DA_1_loading.png)
<img src="/example/DA_1_loading.png" width="600">

```
# Variable loading on discriminant axis 2
axis=2

discriminant_loading(wd = wd, clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var, best.model.number = best.model.number,
                     plot.type = plot.type, plot.width = plot.width, plot.height = plot.height, axis = axis)

```
<img src="/example/DA_2_loading.png" width="600">


```
# Variable loading on discriminant axis 2
axis=3

discriminant_loading(wd = wd, clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var, best.model.number = best.model.number,
                     plot.type = plot.type, plot.width = plot.width, plot.height = plot.height, axis = axis)

```
<img src="/example/DA_3_loading.png" width="600">

<br/>
We can see that variables 1-10 tend to load heavily on at least one discriminant axis, while variables 11-20 do not load heavily on any axes. 
This is expected given that variables 1-10 are simulated to have large differences between group means and small standard deviations,
while variables 11-20 have small differences between group means and large standard deviations.
