# Adjusted predictions in Generalized Estimating Equations (GEEs)

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11296754.svg)](https://doi.org/10.5281/zenodo.11296754) -->

<!-- badges: end -->

This is a Github repository associated with the manuscript "Adjusted Predictions for Generalized Estimating Equations" by [Hui](https://francishui.netlify.app/), [Mueller](https://researchers.mq.edu.au/en/persons/samuel-muller), and [Welsh](https://cbe.anu.edu.au/about/staff-directory/professor-alan-welsh), which is accepted for publicationv in *Biometrics*. As the name suggests, the repository contains code for being able to constructed an adjusted form of prediction for independent cluster GEEs, based on exploiting the assumed working cross-correlations between the current and new observations within the same cluster.

# Getting started

The repository consists of the following set of folders and files:

-   The `simulations` folder contains template scripts `xxx_n25m10.R` for four response types (binary, Gamma, Gaussian, Poisson) to reproduce the simulations in the associated manuscript. The folder also contains a main `simulationfunction.R` script;

-   The `application_spruces` folder contains the `application.R` script to reproducing the application to the sitka spruce growth dataset in the associated manuscript. Note the data are publicly available as part of the [glmtoolbox](https://cran.r-project.org/web/packages/glmtoolbox/index.html) package. **Users interested in an example for how to construct adjusted GEE predictions are recommended to start here**

-   The `R` folder contains `R` scripts for constructing standard and adjusted GEE predictions, associated uncertainty intervals, and `R` scrips used as part of running the simulation study.

    -   Users interested in the "inner workings" of how adjusted GEE predictions are implemented can head to the `predictions.R` function here.

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem, then that is also much appreciated.
3.  Required data files etc...

Alternatively, please contact the corresponding author at [francis.hui\@anu.edu.au](mailto:francis.hui@anu.edu.au)
