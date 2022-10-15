#' Print a summary of the estimated causal mediation model
#'
#' @param fit A \code{stanfit} object obtained from \code{bayesgmed()}
#' @param level "Confidence" level; Defines the limits of the credible intervals.
#' Defaults to .95 (i.e. displays 95\% CIs.)
#' @param pars Parameters to summarize. Defaults to main causal effect estimands
#' similar to the \code{mediation} r package. See Details for more information.
#' @param digits How many decimal points to display in the output. Defaults to 3.
#'
#' @return A \code{data.frame} summarizing the estimated causal
#' mediation model:
#' \describe{
#'  \item{Parameter}{Name of parameters}
#'  \item{Mean}{Mean of parameter's posterior distribution.}
#'  \item{Median}{Median of parameter's posterior distribution.}
#'  \item{SE}{Standard deviation of parameter's posterior distribution.}
#'  \item{ci_lwr}{The lower limit of Credible Intervals.}
#'  \item{ci_upr}{The upper limit of Credible Intervals.}
#'  \item{n_eff}{Number of efficient samples.}
#'  \item{Rhat}{A value of 1.00 suggests model convergency}
#'}
#'
#' @details After estimating a model  with \code{bayesgmed()}, show the estimated results
#' by using \code{bayesgmed_summary(fit)}, where \code{fit} is an object containing
#' the fitted model.
#'
#' The function shows, for each parameter specified with \code{pars},
#' the posterior mean, and limits of the Credible Interval as specified
#' by \code{level}. 
#'
#' \subsection{Parameters}{
#' By default, \code{bayesgmed()} estimates and returns a large number of parameters
#' and their associated standard deviations.
#' However, \code{bayesgmed_summary()} by default only displays a subset of the
#' estimated parameters:
#'
#' \describe{
#'  \item{NDE_control]}{Direct effect estimate when the exposure level is set to the control value.}
#'  \item{NDE_treated]}{Direct effect estimate when the exposure level is set to the treated value.}
#'  \item{NIE_control]}{mediated effect estimate when the exposure level is set to the control value.}
#'  \item{NIE_treated]}{mediated effect estimate when the exposure level is set to the treated value.}
#'  \item{ANDE}{Average direct effect of X on Y}
#'  \item{ANIE}{Average indirect effect of X on Y}
#'  \item{TE}{The total effect of A on Y.}
#'}
#' The user may specify \code{pars = NULL} to display all estimated parameters where all chaines merged.
#' With this argument, \code{bayesgmed_summary()} prints the all parameters defined in the model definations. But the most important ones are:
#'
#' \describe{
#'  \item{alphaZ[]}{Prameter estimate of the confounders (i.e., X -> Y) relationship listed in the order they are specified in covarietes in \code{bayesgmed()}. alpha[1] is the intercept.}
#'  \item{alphaM}{Prameter estimate of the M -> Y relationship.}
#'  \item{alphaA}{Prameter estimate of the A -> Y relationship.}
#'  \item{betaZ}{Prameter estimate of the the confounders (i.e., X -> M) relationship listed in the order they are specified in covarietes \code{bayesgmed()}. beta[1] is the intercept.}
#'  \item{betaA}{Prameter estimate of the A -> M relation ship.}
#'  \item{NDE_control]}{Direct effect estimate when the exposure level is set to the control value.}
#'  \item{NDE_treated]}{Direct effect estimate when the exposure level is set to the treated value.}
#'  \item{NIE_control]}{mediated effect estimate when the exposure level is set to the control value.}
#'  \item{NIE_treated]}{mediated effect estimate when the exposure level is set to the treated value.}
#'  \item{ANDE}{Average direct effect of X on Y}
#'  \item{ANIE}{Average indirect effect of X on Y}
#' #'}

#' To learn more about the additional parameters, refer to the Stan code
#' (\code{cat(get_stancode(fit))}).
#'}
#'
#' @author Belay B. Yimer \email{belaybirlie.yimer@manchester.ac.uk}
#'
#' @export
bayesgmed_summary <- function(
    mod = NULL,
    level = .95,
    pars = c("NDE_control", "NDE_treated", "NIE_control", "NIE_treated", "ANDE", "ANIE","TE"),
    digits = 3
    ){

    # Check that mod is a Stanfit object
    if (!(class(mod) == "stanfit")) stop("Model is not a stanfit object.")

    # Choose which parameters to display
    if (is.null(pars)) pars <- mod@sim$pars_oi  # Return all parameters

    # Obtain model summary from Stanfit
    lower_ci <- .5 - (level/2)
    upper_ci <- .5 + (level/2)

    mod_sum <- rstan::summary(object = mod,
                              probs = c(lower_ci, .5, upper_ci),
                              pars = pars)$summary[,-2]
    # Clean and get post. probs
    if (is.null(dim(mod_sum))) {  # If only one param entered
        mod_sum <- data.frame(t(mod_sum))
        mod_sum <- data.frame(t(apply(mod_sum, 2, round, digits = digits)))
        Names <- pars
    } else {
        mod_sum <- data.frame(mod_sum)
        mod_sum <- data.frame(apply(mod_sum, 2, round, digits = digits))
        Names <- row.names(mod_sum)
        }
    mod_sum$n_eff <- floor(mod_sum$n_eff)
    mod_sum$Parameter <- Names
    mod_sum <- mod_sum[,c(8,1,2,4,3,5,6,7)]
    names(mod_sum) <- c("Parameter", "Mean", "SE", "Median",
                        paste0(lower_ci*100, "%"), paste0(upper_ci*100, "%"),
                        "n_eff", "Rhat")
    row.names(mod_sum) <- NULL
    return(mod_sum)
}