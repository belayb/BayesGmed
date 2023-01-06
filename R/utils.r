#' Print a summary of the estimated causal mediation model
#'
#' @param model A `stanfit` object obtained from `bayesgmed()`.
#' @param level The "confidence" level that defines the limits of the credible intervals (default is .95, i.e. 95% CIs).
#' @param pars The parameters to summarize (defaults to main causal effect estimands similar to the `mediation` package; see Details).
#' @param digits The number of decimal points to display in the output (default is 3).
#'
#' @return A `data.frame` summarizing the estimated causal mediation model, including the following columns:
#'   - `Parameter`: the name of the parameter.
#'   - `Mean`: the mean of the parameter's posterior distribution.
#'   - `Median`: the median of the parameter's posterior distribution.
#'   - `SE`: the standard deviation of the parameter's posterior distribution.
#'   - `ci_lwr`: the lower limit of the credible interval.
#'   - `ci_upr`: the upper limit of the credible interval.
#'   - `n_eff`: the number of efficient samples.
#'   - `Rhat`: a value of 1.00 suggests model convergence.
#'
#' @details After estimating a model with `bayesgmed()`, use `bayesgmed_summary(fit)` to display the estimated results,
#'   where `fit` is an object containing the fitted model. By default, `bayesgmed_summary()` only displays a subset of
#'   the estimated parameters:
#'     - `NDE_control`: direct effect estimate when the exposure level is set to the control value.
#'     - `NDE_treated`: direct effect estimate when the exposure level is set to the treated value.
#'     - `NIE_control`: mediated effect estimate when the exposure level is set to the control value.
#'     - `NIE_treated`: mediated effect estimate when the exposure level is set to the treated value.
#'     - `ANDE`: average direct effect of X on Y.
#'     - `ANIE`: average indirect effect of X on Y.
#'     - `TE`: the total effect of A on Y.
#'   To display all estimated parameters where all chains merged, set `pars = NULL`. This will print all parameters defined
#'   in the model definitions, including the most important ones:
#'     - `alphaZ[]`: parameter estimate of the confounders (i.e., X -> Y) relationship, listed in the order they are specified
#'       in the `covariates` argument of `bayesgmed()`. `alpha[1]` is the intercept.
#'     - `alphaM`: parameter estimate of the M -> Y relationship.
#'     - `alphaA`: parameter estimate of the A -> Y relationship.
#'     - `betaZ`: parameter estimate of the confounders (i.e., X -> M) relationship, listed in the order they are specified
#'       in the `covariates` argument of `bayesgmed()`. `beta[1]` is the intercept 
#'
#' To learn more about the additional parameters, refer to the Stan code
#' (\code{cat(get_stancode(fit))}).
#'}
#'
#' @author Belay B. Yimer \email{belaybirlie.yimer@manchester.ac.uk}
#'
#' @export
bayesgmed_summary <- function(
    model = NULL,
    level = .95,
    pars = c("NDE_control", "NDE_treated", "NIE_control", "NIE_treated", "ANDE", "ANIE","TE"),
    digits = 3
    ){

    # Check that mod is a Stanfit object
    if (!inherits(model, "stanfit")) stop("Model is not a stanfit object.")

    # Choose which parameters to display
    if (is.null(pars)) pars <- model@sim$pars_oi  # Return all parameters

    # Obtain model summary from Stanfit
    lower_ci <- .5 - (level/2)
    upper_ci <- .5 + (level/2)

    model_sum <- rstan::summary(object = model,
                              probs = c(lower_ci, .5, upper_ci),
                              pars = pars)$summary[,-2]
    # Clean and get post. probs
    if (is.null(dim(model_sum))) {  # If only one param entered
        model_sum <- data.frame(t(model_sum))
        model_sum <- data.frame(t(apply(model_sum, 2, round, digits = digits)))
        Names <- pars
    } else {
        model_sum <- data.frame(model_sum)
        model_sum <- data.frame(apply(model_sum, 2, round, digits = digits))
        Names <- row.names(model_sum)
        }
    model_sum$n_eff <- floor(model_sum$n_eff)
    model_sum$Parameter <- Names
    model_sum <- model_sum[,c(8,1,2,4,3,5,6,7)]
    names(model_sum) <- c("Parameter", "Mean", "SE", "Median",
                        paste0(lower_ci*100, "%"), paste0(upper_ci*100, "%"),
                        "n_eff", "Rhat")
    row.names(model_sum) <- NULL
    return(model_sum)
}