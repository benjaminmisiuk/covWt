# covWt
### Spatial covariance weighting in R

This package provides functionality and scripts that were used to perform residual covariance-weighted bagging and validation as presented in the manuscript, *Improved environmental mapping and validation using bagging models with spatially clustered data*. 

Code used to perform the simulations for each validation method presented in the study are found in the "covWt\scripts" directory. Suggest opening the script in R to load in custom data and view code annotations.

## Installation

You can use the `remotes` package to install the package directly from github. Install `remotes` first if you do not have it. 

```
install.packages("remotes")
remotes::install_github("benjaminmisiuk/covWt")
```

## Loading functions
The functions can then be loaded into your R environment.

```
library(covWt)
help(package = 'covWt')
```
