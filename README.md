# BayesGmed - An R package for Bayesian Causal Mediation Analysis using Stan 

# BayesGmed

The R package **BayesGmed** implements parametric mediation analysis using the Bayesian g-formula approch.
In addition to the estimation of causal mediation effects, the package also allows researchers to conduct sensitivity analysis. The methodology behind the R-package and a demonstartion of its application can be found on [arxiv] (https://arxiv.org/abs/2210.08499)
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


# Usage 

The **BayesGmed** R-package currently handles continuous outcome – continuous mediator, binary outcome – binary mediator, continuous outcome – binary mediator, and binary outcome – continuous mediator. 

Suppose we are interested in the causal direct and indirect effect of a single exposure $A$ on a binary outcome $Y$ where we have a single continuous mediator $M$. In addition, assume we have two confounding variables $Z=(Z_1,Z_2)$. The example_data corresponding to this scenerio is included with **BayesGmed**. 

To estimate the direct and indirect of the exposure on the outcome adjusted for confounder, the anlaysis would follow as below. 

``` r
fit1 <- bayesgmed(outcome = "Y", mediator =  "M", treat = "A", covariates = c("Z1", "Z2"), dist.y = "binary",
dist.m = "continuous", link.y = "logit", link.m = "identity", data = example_data)

bayesgmed_summary(fit1)

```

## Contributors

Maintained by Belay Birlie Yimer of the [Centre for Epidemiology Versus Arthritis](https://www.cfe.manchester.ac.uk/), University of Manchester, UK. Co-authors: Mark Lunt, John McBeth, Marcus Beasley, and Gary J Macfarlane. Stan model defination within the package are based on Comment, Leah (2018) Causal inference with the g-formula in Stan. Zenodo.

Pull requests and GitHub issues are welcome.
