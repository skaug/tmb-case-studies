# Case studies in TMB

[TMB](https://github.com/kaskr/adcomp) is software package based on AD (atomatic differentiation) for fitting statistical models. The
present repository contains a collection of worked TMB examples provided for educational purposes. The case studies vary widely in nature and complexity.  


Use
----------------------
* Case studies are  found in separate folders
* Open the `.Rmd` (R Markdown) file in RStudio and push the buttom `knit`

Overview
----------------------

**Model examples**

`simpleExample`
Simple linear regression (your first TMB model). Contains
some basic C++ intro.

`pSplines`
GAM model based on penalized splines


***
**TMB features**

`map_example`
Turning parameters on-and-off in the estimation

`debug_tutorial`
Runtime C++ debugging

***
**Spatial statistics**

`spde`
Spatial modeling with (approximate) Matern covariance functions

`SPDExAR1`
Spatio-temporal modeling

***

**Fisheries**

`SAM`
Stock assessment model

`spatialALK`
Spatial age-length-key model

***

**Wish list**

User contributions are welcome. Topics
of interest are

`Figthing with the C++ compiler`
An example of how to use casting to 
combine `vector`, `matrix` and `array`
objects in various variants. 


