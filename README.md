# MOSS

[![Travis-CI Build Status](https://travis-ci.org/wilsoncai1992/MOSS.svg?branch=master)](https://travis-ci.org/wilsoncai1992/MOSS)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/hagh8vidrdeacr7f?svg=true)](https://ci.appveyor.com/project/wilsoncai1992/MOSS)
[![codecov](https://codecov.io/gh/wilsoncai1992/MOSS/branch/master/graph/badge.svg)](https://codecov.io/gh/wilsoncai1992/MOSS)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MOSS)](http://cran.rstudio.com/web/packages/MOSS/index.html) -->
<!-- [![](http://cranlogs.r-pkg.org/badges/MOSS)](http://cran.rstudio.com/web/packages/MOSS/index.html) [![](http://cranlogs.r-pkg.org/badges/grand-total/MOSS)](http://cran.rstudio.com/web/packages/MOSS/index.html) -->

```R
devtools::install_github('wilsoncai1992/MOSS')
```

## Documentation

* To see all available package documentation:

```R
?MOSS
help(package = 'MOSS')
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

```R
# fit counterfactual survival curve S(t)_A=1
MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e-2, max.iter = 2e2, verbose = FALSE)
MOSS_fit$onestep_curve(g.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                         Delta.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                         ht.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"))
# extract the survival curve S_A=1
MOSS_fit$plot_onestep_curve()
# plot pointwise CI
MOSS_fit$plot_CI_pointwise()
```


## Citation
To cite `MOSS` in publications, please use:
> Cai W, van der Laan MJ (2016). *One-step TMLE for time-to-event outcomes.* Working paper.

## Funding

## Copyright
This software is distributed under the GPL-2 license.

## Community Guidelines
Feedback, bug reports (and fixes!), and feature requests are welcome; file issues or seek support [here](https://github.com/wilsoncai1992/onestep_survival/issues).