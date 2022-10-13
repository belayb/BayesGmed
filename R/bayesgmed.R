#' bayesgmed: Estimate a causal mediation effects
#'
#' Estimates a Bayesian causal mediation model using g-formula approach though Stan.
#' @export
#' @param data A \code{data.frame} or a \code{tibble} object.
#' @param outcome a character string indicating the name of the outcome variable.
#' @param treat a character string indicating the name of the treatment variable. The treatment can be either binary (integer or a two-valued factor) or continuous (numeric).
#' @param mediator a character string indicating the name of the mediator variable.
#' @param covariates a character vector indicating the name of the confounding variables.
#' @param priors A list of named values to be used as the prior scale parameters. See details.
#' @param dist.y a character string indicating the family distribution of the outcome.  E.g., dist.y = "bernoulli" will fit
#' a logistic regression the outcome.
#' @param dist.m a character string indicating the family distribution of the mediator  E.g., dist.m = "bernoulli" will fit
#' a logistic regression the mediator
#' @param link.y a character string indicating the link function to be used for the outcome model.
#' @param link.m a character string indicating the link function to be used for the mediator model.
#' @param control.value The vaue of the treatment variable to be used as a reference level.
#' @param treat.value The vaue of the treatment variable to be used as a treatment level.
#' @param ... Other optional parameters passed to \code{rstan::stan()}.
#'
#' @return An object of S4 class stanfit, with all its available methods.
#'
#' @author Belay Birlie Yimer \email{belaybirlie.yimer@manchester.ac.uk}
#'
#' @details Bayesgmed is the main function for estimating causal mediation effects
#' from several types of data with in the Bayesian framework. We followed the potential
#' outcome framework for effects definition and the package uses the rstan utility functions
#' for exploring the posterior distribution.
#'
#' ## priors
#' Users may pass a list of named values for the priors argument. The values will be used to define
#' the scale parameter of the respective prior distributions. This list may specify some or all of the
#' following parameters:
#' priors <- list(
#'   scale_m = 2.5, scale_y = 2.5,
#'   location_m = 0, location_y = 0,
#'  scale_sd_y = 2.5, scale_sd_m = 2.5)
#'
#'

bayesgmed <- function(outcome, mediator, treat,covariates =NULL,
                     dist.y = "continuous", dist.m ="continuous",
                     link.y ="identity", link.m ="identity",
                     data,  control.value = 0,
                     treat.value = 1, priors = NULL, ...){

  # Check for data
  if (is.null(data)) stop("No data entered")
  if (class(data)[1] == "tbl_df") data <- as.data.frame(data)  # Allow tibbles


  # Check priors
  default_priors <- list(
    scale_m = 2.5, scale_y = 2.5,
    location_m = 0, location_y = 0,
    scale_sd_y = 2.5, scale_sd_m = 2.5
  )
  if (is.null(priors$scale_m)) priors$scale_m <- default_priors$scale_m
  if (is.null(priors$scale_y)) priors$scale_y <- default_priors$scale_y
  if (is.null(priors$location_m)) priors$location_m <- default_priors$location_m
  if (is.null(priors$location_y)) priors$location_y <- default_priors$location_y
  if (dist.y == "continuous"& is.null(priors$scale_sd_y)) priors$scale_sd_y <- default_priors$scale_sd_y
  if (dist.m == "continuous" & is.null(priors$scale_sd_m)) priors$scale_sd_m <- default_priors$scale_sd_m

  # Create a data list for Stan
  stan_data <- list()
  stan_data$X = cbind(1,data[,covariates])
  stan_data$A = data[,treat]
  stan_data$M = data[,mediator]
  stan_data$Y = data[,outcome]
  stan_data$P = ncol(stan_data$X)
  stan_data$N <- length(stan_data$Y)
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
