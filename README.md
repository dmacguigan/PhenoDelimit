# PhenoDelimit
R package to create STRUCTURE-esque barplots and compare species delimitation models with multivariate phenotype data.

This package requires an external program, CLUMPP, which can be [dowloaded here](https://rosenberglab.stanford.edu/clumpp.html). Make sure the CLUMPP executable is in your PATH environmental variable before proceeding. 
For help setting the PATH variable, see these links for [Windows 10](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/) or [Linux/macOS](https://opensource.com/article/17/6/set-path-linux)

To download and install this package
```
install.packages("devtools")
library(devtools)
install_github("dmacguigan/PhenoDelimit")
library(PhenoDelimit)
```

## Example analysis
Data are in "example" directory, so clone the repository first. Data were generated using simulate_data.R in the example folder.
Dataset consists of 20 continuous variables for 300 idividuals divided into four groups.
Data were simulated by drawing from normal disrtibutions.
For each variable, means were allowed to vary between groups, but standard deviations were kept constant.

Variables 1-10 have a large range of group means (between 5 and 30) and small standard deviations (between 1 and 3). These variables are expected to be more informative about group membership.

Variables 11-20 have a small range of group means (between 18 and 22) and large standard deviations (between 3 and 10). Thus, these variables are expected to be less informative about group membership.

![sim_data_boxplots](/example/sim_data_boxplots.png)

<br/>

We will test 6 different delimitation models:
* M1 - model used to simulate the data (4 groups total)
* M2 - merge groups 1 + 2 and 3 + 4 (2 groups total)
* M3 - merge groups 1 + 2 (3 groups total)
* M4 - randomly split group 4 into 2 groups (5 groups total)
* M5 - randomly take samples from groups 3 and 4 and create a new group (5 groups total)
* M6 - keep 4 equal sized groups but randomize individual assignment (4 groups total)

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
```

<br/>

### step 2: summarize CLUMPP
```
wd = "H:/PhenoDelimit/example"
model.numbers = c(1:6)
perc.var = c(70,80,90)

clumpp_results <- read_clumpp_results(wd=wd,
                                      perc.var=perc.var,
                                      model.numbers=model.numbers)
```

<br/>

### step 3: plot H' values to compare delimitation models
```
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
```
![H_plot_example](/example/H_plot_example.png)

Remember, model 1 is the "true" model which generated the data.
So it's comforting to see that model 1 has a near perfect H'.
Models 2-5 are tweaked versions of model 1 with groups merged or split.
Model 6 randomly shuffled the "true" group assignments from model 1.
So no surprise that model 6 has the lowest H'.

<br/>

### step 4: permutation test of significance for H' values - IN DEVELOPMENT
Here we will permute the morphological data and rerun the DAPC and CLUMPP analyses.
We can use this to test the statistical significance of our models.
This may take a while depending on how many models you have and how many permutations you wish to perform.

```
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
```
![example_perm-vs-obs_m1.png](/example/example_perm-vs-obs_m1.png)
![example_perm-vs-obs_m6.png](/example/example_perm-vs-obs_m6.png)

<br/>

In these figure, we can see that the observed H' (red dashed line) for model 1 is signficantly larger than the distribution of permuted H' values (gray bars).
On the other hand, the observed H' (red dashed line) for model 6 is not signficantly different than the permuted H' distribution (gray bars).
This makes sense, since model 1 was used to generated this simulated data, while model 6 was itself a permutation of the underlying data.

![example_obs-minus-perm-mean.png](/example/example_obs-minus-perm-mean.png)

Models 2-5 all have significantly larger H' values than their permuted distributions, but not as large of a differnce as we see for model 1. 
This also makes sense, since models 2-5 are rearrangements of the groups in model 1 used to simulate the data.
But models 2-5 do not reshuffle the entire dataset like model 6.

<br/>


### step 5: bar plots of discriminant analysis assignment probabilities
```
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
![m1_barplot](/example/m1_barplot.png)

<br/>

### step 6: scatter plot or density plot of discriminant axes
```
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
![m1_DA1-DA2_scatter.png](/example/m1_DA1-DA2_scatter.png)

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
![m1_DA1-DA3_scatter.png](/example/m1_DA1-DA3_scatter.png)

<br/>

### step 7: plot discriminant axis loadings and write table of variable contributions and loadings
```
wd = "H:/PhenoDelimit/example/"
clumpp.wd = "H:/PhenoDelimit/example/CLUMPP"
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
<img src="/example/m1_DA_1_loading.png" width="600">

```
# Variable loading on discriminant axis 2
axis=2

discriminant_loading(wd = wd, clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var, best.model.number = best.model.number,
                     plot.type = plot.type, plot.width = plot.width, plot.height = plot.height, axis = axis)

```
<img src="/example/m1_DA_2_loading.png" width="600">


```
# Variable loading on discriminant axis 3
axis=3

discriminant_loading(wd = wd, clumpp.wd = clumpp.wd,
                     best.perc.var = best.perc.var, best.model.number = best.model.number,
                     plot.type = plot.type, plot.width = plot.width, plot.height = plot.height, axis = axis)

```
<img src="/example/m1_DA_3_loading.png" width="600">

<br/>

We can see that variables 1-10 tend to load heavily on at least one discriminant axis, while variables 11-20 do not load heavily on any axes. 
This is expected given that variables 1-10 are simulated to have large differences between group means and small standard deviations,
while variables 11-20 have small differences between group means and large standard deviations.
