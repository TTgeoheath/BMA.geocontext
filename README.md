# About
This package contains four functions corresponding to a two-step approach for addressing residence-based geographic uncertainties. The lacunarity can be
used to determine the upper limits for appropriate buffer selection (gbl).The upper limits is found at the location of maximum curvature (lac.limit).
However, lacunarity assessment only yields the upper-scale limit of the buffer. Even when using buffer sizes within this range, 
variations in the estimated effect sizes in health-environment associations may occur. To address this 
issue, buffer sizes smaller than the upper limit can be considered as a whole, while also taking into account the requirement for exposure-specific buffer 
sizes in multi-exposure models. Multiple-scale Cox models can be developed based on the buffer sizes delineated within the lacunarity-defined limits. 
Bayesian Model Averaging can be extended to Cox models to average effect estimates across all models, thus providing a more robust estimate. Each model is 
weighted by posterior probability. The functions in our package can provide both the coefficients and posterior probabilities for each model (model.pp) 
and an averaged estimate (model.overall). 

# Installation

You can install and load the released version of BMA.geocontext from GitHub with:

```r
library(devtools)
devtools::install_github("TTgeoheath/BMA.geocontext")
library(BMA.geocontext)
```

## Load the test data

We provided two datasets to test our functions. The test_lacunarity can be used to test the function lac.limit for detecting the precise upper limit 
of lacunarity curve. The test_survival can be used to test the model.pp (return the coefficients and posterial coefficient of each model) and 
model.overall ( the averaged results for all models)

```r
# Usage
data(test_lacunarity, package = "BMA.geocontext")
test_lacunarity <- as.data.frame(la)
data(test_survival,  package = "BMA.geocontext")
test_survival <- as.data.frame(test_survial)
# Results
 ---------- test_lacunarity
    box_width lacunarity_values
2          60          2.175266
3          90          1.722891
4         120          1.573135
5         150          1.489790
6         180          1.433522
7         210          1.392854
8         240          1.361289
9         270          1.335002
10        300          1.312525
11        330          1.293059
12        360          1.276228
13        390          1.261002
14        420          1.247146
---------- test_survival
   Status Survival_time Age Income Gender Origin Urbanization  bufND100  bufND200  bufND300  bufND400 
1       0             7  63     44      2      1            4 0.3187178 0.3468641 0.3783087 0.4029252
2       0             6  58     61      2      1            3 0.3615594 0.3767463 0.3918016 0.4149369
3       0             4  42     95      2      1            1 0.3615594 0.3800519 0.3918047 0.4144396
4       0             7  41     83      2      1            2 0.3683210 0.3793667 0.3917254 0.4152211
5       0             8  89     68      1      0            1 0.3683210 0.3812389 0.3917277 0.4154651
6       0             1  80     23      1      1            3 0.3656503 0.3807799 0.3913035 0.4161394
7       0             7  63     83      1      1            1 0.3713914 0.3807370 0.3898436 0.4161140
8       0             4  72     55      2      1            2 0.3685139 0.3805727 0.3880017 0.4162201
9       0             7  96    100      2      1            1 0.3685139 0.3802205 0.3869436 0.4165731
10      0             8  49     90      1      1            2 0.3685139 0.3783150 0.3856445 0.4150494
11      0             3  46     72      2      1            5 0.3761529 0.3783105 0.3879458 0.4144126
12      0             6  53     17      2      1            2 0.3750311 0.3769454 0.3879346 0.4139478
13      0             1  62     10      2      1            1 0.3771406 0.3783995 0.3881606 0.4124978
14      0             2  71     76      1      1            1 0.3615594 0.3758386 0.3885027 0.4144198
15      0             4  62     90      1      1            3 0.3185388 0.3473061 0.3771310 0.4016598
 

```
## Lacunarity analysis
Unsing the funtion gbl(init_boxwidth,end_boxwidth,image, obserwin = Frame(image), raster_type = c("continuous", "binary"));
"init_boxwidth" is the box width for the beginning box put in the corner of a raster surface;
"end_boxwidth" is the box width for your identified last box;
"image" is the image formate of the raster layer;
"obserwin" is the optional observation window. If provided, it will be the intersection of the specified observation window (obserwin) and the 
non-NA pixels in the input image (image);
"raster_type" is the selection of the raster type. It can be either "continuous" or "binary".


```r

# Usage
library(raster);library(maptools);library(RcppRoll)
study_area <- raster(study_area.tif)
im_area <- as.im.Raster(study_area)
LA <- gbl(30,10000,im_area,raster_type="continuous")

		
# Results
> LA
Function value object (class ‘fv’)
for the function BoxWidth -> expression(lacunarity_values[gb])
..............................................................................
                  Math.label      Description                                 
box_width         BoxWidth        the box width of gliding box                
lacunarity_values Lacunarityvalue Lacuanrity vaule corresponding to box widths
..............................................................................
Default plot formula:  .~.x
where “.” stands for ‘lacunarity_values’
Recommended range of argument box_width: not specified
Available range of argument box_width: [60, 10260]
## Critical soil moisture and maximum bulk density (Proctor test)
> View(LA)
   box_width lacunarity_values
2          60          2.175266
3          90          1.722891
4         120          1.573135
5         150          1.489790
6         180          1.433522
7         210          1.392854
8         240          1.361289
9         270          1.335002
10        300          1.312525
11        330          1.293059
12        360          1.276228
13        390          1.261002
14        420          1.247146


```
## Compute the flattened points (upper limit)
The function lac.limit(data, box_sizes, lacunarity_values,start = list(a = 1, b = -1, c = 1)) computes the precise point where the lacunarity curve
flattened out.
"data" is a dataframe includes the boxsizes of the gliding boxes and their corresponding lacuanrity value;
"box_sizes" The increased box widths of gliding boxes in lacunarity analysis;
"lacunarity_values" is the corresponding lacunarity values of the increased box sizes given by lacunarity analysis;
"start" is the start value of the curvature simulation. This function used a power function a * x^b + c to simulate the lacunarity curvature. Through
simulation, the flattened point can be obtained by computing the maximum curvature. A start = list(a = 1, b = -1, c = 1) was set as default. Users can 
adjusted it based on their data.
```
# Usage
  maxcurv <- lac.limit(test_lacunarity, "box_width", "lacunarity_values")
# Results
------------------------------------
        Maximum curvature point 
 
f(x) = a * x^b + c 
critical x:  984.3049 
critical y:  1.118582
method: perpendicular distances 
--------------------------------------

```

## Bayesian model averaging pooled cox regression

add later

```r


# Usage

# Results
add later
 
# Citation and references

add later
