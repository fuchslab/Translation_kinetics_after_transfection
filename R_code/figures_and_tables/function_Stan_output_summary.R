#' constructs a summary matrix for one stanfit object
#'
#' @param stanfit_object stanfit object as return by the function rstan::stan()
#' @param param the parameter names as contained in the stanfit_object
#' @param param_names the parameter names as to be included as row names of the 
#' summary matrix for printing
#' @param true_values the true parameter values if available, default: NA
#'
#' @return the summary matrix

tab_Stan_output_summary <- function(stanfit_object, param, param_names, 
                                    true_values = NA){
  # the summary.stanfit method returns the following columns:
  # "mean", "se_mean", "sd", "2.5%", "50%", "97.5%",  "n_eff", "Rhat" 
  stan_summary <- summary(stanfit_object,
                          pars = param,
                          probs = c(0.025, 0.5, 0.975))$summary
  # extract the MCMC samples from the stanfit objects and convert them to a matrix
  samples <- do.call(cbind, rstan::extract(stanfit_object, pars = param))
  # calculate the coefficient of variation of the posterior sample for each paramter
  cv_vec <- apply(samples, 2, sd) / apply(samples, 2, mean) 
  # contruct the final summary matrix 
  if(all(is.na(true_values))){
    stan_summary <- cbind("mean" = stan_summary[ , "mean"], 
                          "c.v." = cv_vec,
                          stan_summary[ , 4:8])  
    # adjust column and row names for printing
    colnames(stan_summary)[3:7] <-  c("$2.5\\%$", "$50\\%$", 
                                      "$97.5\\%$", "$n_\\text{eff}$", "$\\hat{R}$")
  }else{
    stan_summary <- cbind("true value" = true_values, 
                          "mean" = stan_summary[ , "mean"], 
                          "c.v." = cv_vec,
                          stan_summary[ , 4:8])  # adjust column and row names for printing
    colnames(stan_summary)[4:8] <-  c("$2.5\\%$", "$50\\%$", 
                                      "$97.5\\%$", "$n_\\text{eff}$", "$\\hat{R}$")
  }

  rownames(stan_summary) <- param_names
  stan_summary
}