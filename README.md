# onestep.survival

<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/onestep.survival)](http://cran.rstudio.com/web/packages/onestep.survival/index.html) -->
<!-- [![](http://cranlogs.r-pkg.org/badges/onestep.survival)](http://cran.rstudio.com/web/packages/onestep.survival/index.html) [![](http://cranlogs.r-pkg.org/badges/grand-total/onestep.survival)](http://cran.rstudio.com/web/packages/onestep.survival/index.html) -->
<!-- [![Travis-CI Build Status](https://travis-ci.org/wilsoncai1992/onestep.survival.svg?branch=master)](https://travis-ci.org/wilsoncai1992/onestep.survival) -->

<img style="float: left;margin:0 5rem 0 0" src="http://media.web.britannica.com/eb-media/29/76829-050-CD9C4B43.jpg" width="30%" height="30%">
<br>
<!-- <img style="float: left;margin:0 5rem 0 0" src="http://www.feenixx.com/space-exploration/posters/First_Step_on_Moon_Poster.jpg" width="30%" height="30%">
<br>
 -->
The `onestep.survival` R package is a tool for estimating counterfactual survival curve under static or dynamic interventions on treatment (exposure), while at the same time adjust for *measured* counfounding. Targeted Maximum Likelihood Estimate (TMLE) approach is employed to create a doubly robust and semi-parametrically efficient estimator. Machine Learning algorithms (SuperLearner) are implemented to all stages of the estimation.

Currently implemented **estimators** include:

1. One-step TMLE for the whole survival curve
2. One-step TMLE for survival at a specific end point
3. Iterative TMLE for survival at a specific end point

## Installation

To install the development version (requires the `devtools` package):

```R
install.packages('SuperLearner')
install.packages('tmle')
install.packages('Matrix')
devtools::install_github('wilsoncai1992/survtmle')
devtools::install_github('wilsoncai1992/onestep.survival')
```

## Documentation

* To see all available package documentation:

```R
?onestep.survival
help(package = 'onestep.survival')
```

## Brief overview

### Data structure

The data input of all methods in the package should be an `R` `data.frame` in the following survival long data format:

```R
#   ID W A T.tilde delta
# 1  1 0 0      95     1
# 2  2 1 1       1     0
# 3  3 0 0     215     1
# 4  4 1 1      15     1
# 5  5 0 0      73     1
# 6  6 0 0      15     1
```

## Citation
To cite `onestep.survival` in publications, please use:
> Cai W, van der Laan MJ (2016). *One-step TMLE for time-to-event outcomes.* Working paper.

## Funding

## Copyright
This software is distributed under the GPL-2 license.

## Community Guidelines
Feedback, bug reports (and fixes!), and feature requests are welcome; file issues or seek support [here](https://github.com/wilsoncai1992/onestep_survival/issues).