# BayesGmed - An R package for Bayesian Causal Mediation Analysis using Stan 

# BayesGmed

The R package **BayesGmed** implements parametric mediation analysis using the Bayesian g-formula approch.
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
devtools::install_github("belayb/BayesGmed")
```

## Contributors

Maintained by Belay Birlie Yimer of the [Centre for Epidemiology Versus Arthritis](https://www.cfe.manchester.ac.uk/), University of Manchester, UK. Co-authors: Mark Lunt, John McBeth, Marcus Beasley, and Gary J Macfarlane. Stan model defination within the package are based on (https://doi.org/10.5281/zenodo.1285276).

Pull requests and GitHub issues are welcome.
