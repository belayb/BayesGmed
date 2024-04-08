# BayesGmed - An R package for Bayesian Causal Mediation Analysis using Stan 

# BayesGmed

The R package **BayesGmed** implements parametric mediation analysis using the Bayesian g-formula approach.
In addition to the estimation of causal mediation effects, the package also allows researchers to conduct sensitivity analysis. The methodology behind the R-package and a demonstration of its application can be found on [arxiv] (https://arxiv.org/abs/2210.08499)
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

Suppose we are interested in the causal direct and indirect effect of a single exposure $A$ on a binary outcome $Y$ where we have a single binary mediator $M$. In addition, assume we have two confounding variables $Z=(Z_1,Z_2)$. The example_data corresponding to this scenario is included with **BayesGmed**. 

To estimate the direct and indirect of the exposure on the outcome adjusted for confounder, the anlaysis would follow as below. 

``` r
data(example_data)
fit1 <- bayesgmed(outcome = "Y", mediator =  "M", treat = "A", covariates = c("Z1", "Z2"), dist.y = "binary",
dist.m = "binary", link.y = "logit", link.m = "logit", data = example_data)

bayesgmed_summary(fit1)

```

The above model assumes the following structure for the outcome and mediator model 

\eqn{logit( P(Y_{i} = 1|A_{i},M_{i},\mathbf{Z}_{i}) ) = \alpha_{0} + \mathbf{\alpha}_{Z}^{'}\mathbf{Z}_{i} + \alpha_{A}A_{i} + \alpha_{M}M_{i},}

\eqn{logit(P(M_{i} = 1| A_{i},\mathbf{Z}_{i} ) ) = \beta_{0} + \mathbf{\beta}_{Z}^{'}\mathbf{Z}_{i} + \beta_{A}A_{i}.}

Note: BayesGmed currently does not handle interaction. 
### Priors

Currently, a multi-normal, $MVN(\text{location}, \text{scale})$, prior is assigned to all regression parameters where the location and scale parameters are fixed to the following default values. The user can change the location and scale parameters by passing the location and scale parameters of the priors as a list as below 

``` r
priors <- list(scale_m = 2.5*diag(P+1), 
               scale_y = 2.5*diag(P+2),
               location_m = rep(0, P+1), 
               location_y = rep(0, P+2),
               scale_sd_y = 2.5, 
               scale_sd_m = 2.5)
```
where $P$ is the number of covariates (including the intercept) in the mediator/outcome model but excluding the exposure and mediator. For the residual standard deviation, a half normal prior distribution with mean zero and scale parameter scale_sd_* is assumed. Currently, the user can only change the scale and location parameters. 

## References 
-  McCandless, L.C. and J.M. Somers, Bayesian sensitivity analysis for unmeasured confounding in causal mediation analysis. Statistical Methods in Medical Research, 2019. 28(2): p. 515-531.
- Comment, L., Coull, B. A., Zigler, C., & Valeri, L. (2019). Bayesian data fusion for unmeasured confounding. arXiv preprint arXiv:1902.10613.
- Imai, K., L. Keele, and D. Tingley, A general approach to
causal mediation analysis. Psychological methods, 2010. 15(4):
- Keil, A.P., et al., A Bayesian approach to the g-formula.
Statistical methods in medical research, 2018. 27(10): p.
3183-3204.



## Contributors

Maintained by Belay Birlie Yimer of the [Centre for Epidemiology Versus Arthritis](https://www.cfe.manchester.ac.uk/), University of Manchester, UK. Co-authors: Mark Lunt, John McBeth, Marcus Beasley, and Gary J Macfarlane. Stan model definition within the package are based on Comment, Leah (2018) Causal inference with the g-formula in Stan.

## Note

The package is under development. Hence, pull requests and GitHub issues are welcome. Any use of the package has to be done wth care. 

<!-- badges: start -->
  [![R-CMD-check](https://github.com/belayb/BayesGmed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/belayb/BayesGmed/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
