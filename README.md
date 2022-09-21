# BayesGmed - An R package for Bayesian Causal Mediation Analysis using Stan 

# BayesGmed

The R package **BayesGmed** implements parametric mediation analysis using the Bayesian g-formula approch (based on 10.5281/zenodo.1285275, 10.1177/0962280217729844).
In addition to the estimation of causal mediation effects, the package also allows researchers to conduct sensitivity analysis.

The package’s source code is hosted on [GitHub](https://github.com/belayb/BayesGmed/). More information can be found on the **BayesGmed**’s Vignette .

# Install

## Requirements

Please ensure you have the latest version of R installed. Windows users
may need to install RTools (more information on the [RStan
website](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows)),
OS X users may need to install XCode


## Install from GitHub

**BayesGmed** is developed on GitHub, and users may obtain the latest (development) version from GitHub directly.

The latest development version of **BayesGmed** requires
[devtools](https://cran.r-project.org/package=devtools) for
installation. If you don’t have the devtools package installed in R,
first run this line:

``` r
install.packages("devtools")
```

Then proceed to install **BayesGmed** from GitHub:

``` r
devtools::install_github("belayb/BayesGmed", args = "--preclean")
```

# Information

Please contact the author of the package for questions and suggestions.
I recommend creating a new issue on GitHub.

# Citation

If you use this software, please cite it:
