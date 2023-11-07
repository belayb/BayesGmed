#' Estimate a causal mediation effects
#'
#' Estimates various quantities for causal mediation analysis using 'Stan'.
#' @export
#' @param data A \code{data.frame} or a \code{tibble} object.
#' @param outcome a character string indicating the name of the outcome variable.
#' @param treat a character string indicating the name of the treatment variable. The treatment variable is considered binary and should be coded as 0 for control and 1 for treated.
#' @param mediator a character string indicating the name of the mediator variable.
#' @param covariates a character vector indicating the name of the confounding variables.
#' @param priors A list of named values for the prior scale parameters. See details.
#' @param dist.y a character string indicating the family distribution of the outcome.  E.g., dist.y = "bernoulli" will fit
#' a logistic regression for the outcome.
#' @param dist.m a character string indicating the family distribution of the mediator, E.g., dist.m = "bernoulli" will fit
#' a logistic regression model for the mediator
#' @param link.y a character string indicating the link function to be used for the outcome model.
#' @param link.m a character string indicating the link function to be used for the mediator model.
#' @param ... Other optional parameters passed to \code{rstan::stan()}.
#'
#' @return An object of 'S4' class 'stanfit',  with all its available methods.
#'
#' @author Belay Birlie Yimer \email{belaybirlie.yimer@manchester.ac.uk}
#'
#' @details This is the main function for estimating causal mediation effects
#' from several types of data within the Bayesian framework. We followed the potential
#' outcome framework for effects definition and the package uses the 'rstan' utility functions
#' for exploring the posterior distribution.
#'
#' ## priors
#' Users may pass a list of named values for the priors argument. The values will be used to define
#' the scale parameter of the respective prior distributions. This list may specify some or all of the
#' following parameters:
#' priors <- list(
#'   scale_m = 2.5*diag(P_m) scale_y = 2.5*diag(P_y),
#'   location_m = rep(0, P_m) location_y = rep(0, P_y),
#'  scale_sd_y = 2.5, scale_sd_m = 2.5)
#' where P_m is the number of regression parameters (including the intercept) in the mediator model and
#' P_y is the number of regression parameters in the outcome model.    
#'
#' @references 
#' 1. McCandless, L.C. and J.M. Somers, \emph{Bayesian sensitivity analysis for unmeasured confounding in causal mediation analysis.} Statistical Methods in Medical Research, 2019. (28)(2): p. 515-531.
#' 2. Comment, L., Coull, B. A., Zigler, C., and Valeri, L. (2019). Bayesian data fusion for unmeasured confounding. arXiv preprint arXiv:1902.10613.
#' 
#' @examples
#' \dontrun{
#' ## Run example using the example_data
#' data(example_data)
#' fit1 <- bayesgmed(outcome = "Y", mediator =  "M", treat = "A", covariates = c("Z1", "Z2"),
#'  dist.y = "binary", dist.m = "binary", link.y = "logit", link.m = "logit", data = example_data)
#' bayesgmed_summary(fit1)
#' 
#' # With priors
#' priors <- list(scale_m = 2.5*diag(P+1), 
#               scale_y = 2.5*diag(P+2),
#               location_m = rep(0, P+1), 
#               location_y = rep(0, P+2))
#' fit1 <- bayesgmed(outcome = "Y", mediator =  "M", treat = "A", covariates = c("Z1", "Z2"), 
#' dist.y = "binary", dist.m = "binary", link.y = "logit", link.m = "logit", priors = priors,  
#' data = example_data)
#' bayesgmed_summary(fit1)
#' }

bayesgmed <- function(outcome, mediator, treat,covariates =NULL,
                     dist.y = "continuous", dist.m ="continuous",
                     link.y ="identity", link.m ="identity",
                     data,  priors = NULL, ...){

  # Check for data
  if (is.null(data)) stop("No data entered")
  if (class(data)[1] == "tbl_df") data <- as.data.frame(data)  # Allow tibbles
  # Create a data list for Stan
  stan_data <- list()
  stan_data$X = cbind(1,data[,covariates])
  stan_data$A = data[,treat]
  stan_data$M = data[,mediator]
  stan_data$Y = data[,outcome]
  stan_data$P = ncol(stan_data$X)
  stan_data$N <- length(stan_data$Y)

# Check priors
  default_priors <- list(
    scale_m = 2.5*diag(stan_data$P + 1), scale_y = 2.5*diag(stan_data$P + 2),
    location_m = rep(0, stan_data$P + 1), location_y = rep(0, stan_data$P + 2),
    scale_sd_y = 2.5, scale_sd_m = 2.5
  )
  if (is.null(priors$scale_m)) priors$scale_m <- default_priors$scale_m
    if (!is.null(priors$scale_m)&dim(priors$scale_m)[1] != stan_data$P + 1) stop("Not all priors supplied for the mediator model")
  if (is.null(priors$scale_y)) priors$scale_y <- default_priors$scale_y
      if (!is.null(priors$scale_y)&dim(priors$scale_y)[1] != stan_data$P + 2) stop("Not all priors supplied for the outcome model")
  if (is.null(priors$location_m)) priors$location_m <- default_priors$location_m
        if (!is.null(priors$scale_y)&length(priors$location_m) != stan_data$P + 1) stop("Not all priors supplied for the mediator model")
  if (is.null(priors$location_y)) priors$location_y <- default_priors$location_y
          if (!is.null(priors$scale_y)&length(priors$location_y) != stan_data$P + 2) stop("Not all priors supplied for the outcome model")
  if (dist.y == "continuous"& is.null(priors$scale_sd_y)) priors$scale_sd_y <- default_priors$scale_sd_y
  if (dist.m == "continuous" & is.null(priors$scale_sd_m)) priors$scale_sd_m <- default_priors$scale_sd_m

  stan_data <- append(stan_data, priors)


  
  if (dist.y == "continuous" & dist.m == "continuous"){
    out <- rstan::sampling(stanmodels$NY_NM_single, data = stan_data, ...)
  }
  else if (dist.y == "binary" & dist.m == "continuous"){
    out <- rstan::sampling(stanmodels$BY_NM_single, data = stan_data, ...)
  }
  else if (dist.y == "binary" & dist.m == "binary"){
    out <- rstan::sampling(stanmodels$BY_BM_single, data = stan_data, ...)
  }
  else if (dist.y == "continuous" & dist.m == "binary"){
    out <- rstan::sampling(stanmodels$NY_BM_single, data = stan_data, ...)
  }
  else {
    stop("model not supported")
  }

  return(out)

}
