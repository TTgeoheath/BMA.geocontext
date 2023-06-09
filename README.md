# About
This R package contains four functions corresponding to a two-step approach for addressing residence-based geographic uncertainties using buffer analyses.
The lacunarity is used to determine the upper buffer limits (gbl()). The upper limit is determined through maximizing the curvature of the lacunarity 
curve (lac.limit()). However, lacunarity assessment only yields the upper-scale limit of the buffer. Even when using buffer sizes within this range, 
variations in the estimated effect sizes in the associations may occur. To address this issue, buffer sizes smaller than the upper limits are used within 
a Bayesian model averaging framework. Multi-scale Cox models are developed based on the buffer sizes within the lacunarity-defined upper limits. We 
extended Bayesian Model Averaging to the Cox regression model to average effect estimates across all multi-scale models, thus providing a more robust 
estimate. Each model is weighted by its posterior probability. The functions provide the coefficients and posterior probabilities for each model
(model.pp()) and an averaged estimate (model.overall()). 

# Installation

You can install and load the released version of BMA.geocontext from GitHub with:

```r
library(devtools)
devtools::install_github("TTgeoheath/BMA.geocontext")
library(BMA.geocontext)
```

## Load the test data
We provided two artificial test datasets. The data "test_lacunarity" represents gridded exposure surfaces and is used to test the function lac.limit() for 
detecting the precise upper limit of the lacunarity curve. The data “test_survival” is used to test the model.pp() and returns the coefficients and 
posterior coefficient of each model. The function model.overall() provides the averaged results across  the models. 

The variables in "test_survival" are as follows: 1) Survival_time: the death,survival, cencoring time in years of the participants. 2) Status: the 
participates expericed death (1) or not (0) during the period. 3)The age of participates at the baseline. 4) Income: the yearly household income (Euro in 
k). 5) Gender: the gender of participants at the basedline (male=1, felmale=2). 6) Origin: the origin of perticipants (Dutch=1, non-Dutch=0). 
7)Urbanization: the urbanization levels of the municipality where the participants live in. 8) bufDN100...900: The NDVI was delineated by buffer sizes 
from 100 to 900 meters. 9) bufPM100...400: The PM2.5 was delineated by buffer sizes from 100 to 400 meters. 10) bufno100...900: The noise was delineated 
by buffer sizes from 100 to 900 meters.

```r
# Usage
data(test_lacunarity, package = "BMA.geocontext")
test_lacunarity <- as.data.frame(la)
data(test_survival,  package = "BMA.geocontext")
test_survival <- as.data.frame(test_survial)
# Results

> head(test_lacunarity)
 ---------- test_lacunarity---------------------------
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


>head(test_survival)
--------------------------------------test_survival(sample)------------------------------------------
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
---------------------------------------------------------------------------------------------------------

```
## Lacunarity analysis
The function gbl(init_boxwidth,end_boxwidth,image, obserwin = Frame(image), raster_type = c("continuous", "binary")) calculates the lacunarity values against increased box width.

The parameter "init_boxwidth" refers to the box width, while the parameter “end_boxwidth” refers to the box width for your identified last box; "image" is the the raster layer input; "obserwin" is the optional observation window. If provided, it will be the intersection of the specified observation window (obserwin) and the non-NA pixels in the input image (image); "raster_type" is the selection of the raster type. It can be either "continuous" or "binary".

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

The parameter "data" is a dataframe includes the boxsizes of the gliding boxes and their corresponding lacuanrity value;
"box_sizes"  is the increased box widths of gliding boxes in lacunarity analysis;
"lacunarity_values" is the corresponding lacunarity values of the increased box sizes given by lacunarity analysis;
"start" is the start value of the curvature simulation. This function used a power function a * x^b + c to simulate the lacunarity curvature. Through
simulation, the flattened point can be obtained by computing the maximum curvature. A start = list(a = 1, b = -1, c = 1) was set as default. Users can 
adjusted it based on their data.
```
# Usage
  box_width <- test_lacunarity[,c("box_width"]
  lacunarity_values <- test_lacunarity[,c("lacunarity_values"]
  maxcurv <- lac.limit(test_lacunarity, box_width, lacunarity_values)
# Results
> print(maxcurv)
------------------------------------
        Maximum curvature point 
 
f(x) = a * x^b + c 
critical x:  984.3049 
critical y:  1.118582
method: perpendicular distances 
--------------------------------------

```

## Bayesian model averaging pooled cox regression
The function model.pp(data,surv_time,status,covariates,exposure_colnames) calculates the coefficients (i.e., beta, p-value, z-value, standard error, 
hazard ratio, confidence interval), posterior possibility, and BIC for each cox model. The function 
model.overall(data,surv_time,status,covariates,exposure_colnames) computes the averaged results of all models.

The parameter "data" is a data.frame containing the survival time, status, covariates, and exposures delineated by different geographic contexts; 
"surv_time" is the column name of survival time; "status" is the column name of status; "covariates" is a character vector of column names for the fixed 
covariates; "exposure_colnames" is a character vector of column names for the exposure variables; the column names of exposures should be the combination 
of prefixes and buffer sizes, and the same exposures should have the same prefixes (e.g., the column names of green space for buffer sizes 100, 200, 300, 
400 should be  bufnd100, bufnd200, bufnd300, and bufnd400). Note, the prefixes can only be letters, e.g., bufnd, bufpm, bufno.

```r
# Usage
library(survival)
cova <- colnames(test_survival[,c(3:7)])
exposure_col<- colnames(test_survival[,c(8:29)])
pp <- model.pp(test_survival,"Survival_time","Status",cova,exposure_col)
overall <- model.overall(test_survival,"Survival_time","Status",cova,exposure_col)
# Results
> head (pp)
-----------------------------------each model (sample)----------------------------------------------
               V1       V2       V3           bufND_Beta            bufPM_Beta            bufno_Beta
model.1  bufND100 bufPM100 bufno100 -0.00888266705382276  -0.00182208209835566 -2.83349355781883e-05
model.2  bufND200 bufPM100 bufno100   0.0258314756399133 -0.000828747328277729  7.25383170306356e-05
model.3  bufND300 bufPM100 bufno100   0.0253558862553519 -0.000768475213503715  6.65892186767771e-05
model.4  bufND400 bufPM100 bufno100   0.0236506924168014 -0.000767544763627434  6.00001470833901e-05
model.5  bufND500 bufPM100 bufno100   0.0512362723196588   0.00030141051876061  0.000128084837236574
model.6  bufND600 bufPM100 bufno100       0.063984399343  0.000877028981516327  0.000160007240942227
                  bufND_HR          bufPM_HR          bufno_HR          bufND_SE           bufPM_SE
model.1  0.991156667282409 0.998179576885476 0.999971665465852 0.105783392938941 0.0104086719088843
model.2   1.02616799959533 0.999171595987942    1.000072540948 0.116731710865855  0.010624773100112
model.3   1.02568008102159  0.99923181998795  1.00006659143579 0.125497289689191 0.0108454692523499
model.4   1.02393258799771 0.999232749723506  1.00006000194713 0.132955491839692 0.0110686501415377
model.5    1.0525715574293  1.00030145594748  1.00012809304045 0.140284161662206 0.0113096282139341
model.6   1.06607576714939  1.00087741368389  1.00016002004278 0.148279522803943  0.011566840950248
                    bufno_SE             bufND_z               bufPM_z               bufno_z
model.1  0.00204524327580101 -0.0839703360521808    -0.175054234998071   -0.0138540661218363
model.2   0.0020441679689513   0.221289274767831   -0.0780014142861075    0.0354854973428868
model.3  0.00204451809228398    0.20204329765327   -0.0708567970295252    0.0325696402140363
model.4  0.00204583541195505   0.177884283601595   -0.0693440260386445    0.0293279443364667
model.5  0.00204776068059547   0.365232052660598    0.0266507893150063    0.0625487335753162
model.6  0.00205119483064046    0.43151203978179    0.0758226887780906    0.0780068468153593

>View(overall)
--------------------------------------overall------------------------------------------------------
                               [,1]
bufND_Beta              3.014981e-02
bufPM_Beta             -4.498268e-04
bufno_Beta             -6.455892e-05
bufND_HR                1.030918e+00
bufPM_HR                9.995515e-01
bufno_HR                9.999356e-01
bufND_SE                1.462310e-01
bufPM_SE                1.224808e-02
bufno_SE                2.310743e-03
bufND_z                 1.942189e-01
bufPM_z                -3.641898e-02
bufno_z                -2.672774e-02
bufND_P_value           8.261673e-01
bufPM_P_value           9.154665e-01
bufno_P_value           8.318294e-01
bufND_lowerci_overall   7.744102e-01
bufPM_lowerci_overall   9.758432e-01
bufno_lowerci_overall   9.954172e-01
bufND_higherci_overall  1.375459e+00
bufPM_higherci_overall  1.023838e+00
bufno_higherci_overall  1.004475e+00
------------------------------------end--------------------------------------------------------------
 
# Citation and references
Labib, S. M., Lindley, S. & Huck, J. J. Scale effects in remotely sensed greenspace metrics and how to mitigate them for environmental health exposure
assessment

Mandelbrot, B. B. A Fractal’s Lacunarity, and how it can be Tuned and Measured. in Fractals in Biology and Medicine (1994)

Plotnick, R. E., Gardner, R. H., Hargrove, W. W., Prestegaard, K. & Perlmutter, M. Lacunarity analysis: A general technique for the analysis of spatial 
patterns. Phys Rev E 53, 5461–5468 (1996)

Dong, P. Test of a new lacunarity estimation method for image texture analysis. Int J Remote Sens 21, 3369–3373 (2000)

Volinsky, C. T., Madigan, D., Raftery, A. E. & Kronmal, R. A. Bayesian Model Averaging in Proportional Hazard Models: Assessing the Risk of a Stroke.
RStat Soc Ser C Appl Stat 46, 433–448 (1997)
