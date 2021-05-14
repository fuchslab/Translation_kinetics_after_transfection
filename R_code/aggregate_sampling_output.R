# This R script aggregates the Stan sampling results for all observed 
# trajectories separately for each dataset and model type.
library(rstan)

input_args <- commandArgs(trailingOnly = TRUE)

# check input arguments --------------------------------------------------------
if(input_args[1] %in% c("experimental_data_eGFP", "experimental_data_d2eGFP",
                        "simulated_data_dataset1_with_error",
                        "simulated_data_dataset2_with_error",
                        "simulated_data_dataset1_no_error",
                        "simulated_data_dataset2_no_error")){
  dataset <- input_args[1]
} else {
  stop(paste0("Unknown dataset '",  input_args[1],
              "' supplied. Eligible first arguments are 
              [experimental_data_eGFP, experimental_data_d2eGFP, 
              simulated_data_dataset1_with_error, simulated_data_dataset2_with_error, 
              simulated_data_dataset1_no_error, simulated_data_dataset2_no_error]."),
       sep = "")
}


if(input_args[2] %in% c("SDE", "ODE")){
  model_type <- input_args[2]
} else {
  stop(paste0("Unknown model type '",  input_args[1],
              "' supplied. Eligible second arguments are [SDE, ODE]."),
       sep = "")
}


# pre-allocate memory -----------------------------------------------------------
N <- 100
# define matrix to store diagnostics for each stanfit object
names_diagnostics_overall <- c("num_div_trans", "num_treedepth_exc", 
                               "min_BFMI", "mean_BFMI", "num_BFMI_below_02", 
                               "max_total_time", "mean_total_time")
mat_diagnostics_overall <- matrix(data = NA, nrow = N,
                                  ncol = length(names_diagnostics_overall))
colnames(mat_diagnostics_overall) <- names_diagnostics_overall

# split input argument `dataset` at underscores
dataset_split <- strsplit(dataset, '_')[[1]]

if(model_type == "SDE"){
  if(dataset_split[1] == "simulated" & dataset_split[4] == "no"){# no measurement error
    names_param <- c("theta[1]", "theta[2]", "theta[3]", "m0", "scale", 
                     "prod_theta2_m0", "prod_theta2_scale", 
                     "prod_m0_scale", "prod_theta2_m0_scale", "x[180,1]")
  }else{
    names_param <- c("theta[1]", "theta[2]", "theta[3]", "m0", "sigma", "scale", 
                     "offset", "prod_theta2_m0", "prod_theta2_scale", 
                     "prod_m0_scale", "prod_theta2_m0_scale", "x[180,1]", "x[180,2]")
  }
}else if(model_type == "ODE"){
  names_param <- c("theta[1]", "theta[2]", "theta[3]", "m0", "sigma", "scale", 
                   "offset", "t0", "prod_theta2_m0", "prod_theta2_scale", 
                   "prod_m0_scale", "prod_theta2_m0_scale", "x2_sim[180]")
  
}

# define array to store the summary statistics specific for each parameter of each stanfit object
colnames_summary <- c("mean", "se_mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat",
                      "hpdi_l", "hpdi_u")
arr_summary <- array(data = NA, 
                     dim = c(N, length(names_param), length(colnames_summary)),
                     dimnames = list(obs_index = 1:N, parameters = names_param,
                                     summary_statistic = colnames_summary))

# calculate summary statistics for each stanfit object -------------------------
for (obs_index in 1:N){ # for each trajectory
  path_stanfit_object <- paste0("intermediate_output_files/stanfit_objects/",
                                dataset , "_", model_type, "/stanfit_object_",
                                obs_index, ".rds", sep = "")
  stanfit_object <- readRDS(path_stanfit_object)
  
  mat_diagnostics_overall[obs_index, "num_div_trans"] <- 
    get_num_divergent(stanfit_object)
  mat_diagnostics_overall[obs_index, "num_treedepth_exc"] <- 
    get_num_max_treedepth(stanfit_object)
  bmfi <- get_bfmi(stanfit_object)
  mat_diagnostics_overall[obs_index, "min_BFMI"]  <- min(bmfi)
  mat_diagnostics_overall[obs_index, "mean_BFMI"] <- mean(bmfi)
  mat_diagnostics_overall[obs_index, "num_BFMI_below_02"]  <- sum(bmfi < 0.2)
  total_time_elapsed <- rowSums(get_elapsed_time(stanfit_object))
  mat_diagnostics_overall[obs_index, "max_total_time"]  <- max(total_time_elapsed)
  mat_diagnostics_overall[obs_index, "mean_total_time"]  <- mean(total_time_elapsed)
  
  s <- summary(stanfit_object, pars = names_param, 
               probs = c(0.025, 0.5, 0.975))$summary
  s <- round(s, digits = 3)
  s[,"n_eff"] <- round(s[,"n_eff"], digits = 0)
  arr_summary[obs_index, , 
              c("mean", "se_mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")] <- s
  # calculate the highest probability density intervals (HPDIs)
  # extract sample from stanfit object
  samples <- as.matrix(stanfit_object, pars = c(names_param, "lp__"))
  # order sample according to value of the log posterior
  ordered_samples <- samples[order(samples[,"lp__"]),]
  # determine the fraction of the sample with the (95%) highest log posterior
  # and exclude the column of the log posterior
  hpd_samples <- ordered_samples[-(1:(floor(dim(ordered_samples)[1]*0.05))),
                                 !dimnames(ordered_samples)$parameters %in% c("lp__")]
  # determine the range of values with the highest log posterior for each parameter
  arr_summary[obs_index, , c( "hpdi_l", "hpdi_u")] <- 
    round(t(apply(hpd_samples, 2, range)), digits = 3)
}


# save aggregated output -------------------------------------------------------

path_aggregated_output_folder  <- 
  paste0("intermediate_output_files/aggregated_output/", dataset, sep = "")
dir.create(file.path(getwd(), path_aggregated_output_folder), 
           showWarnings = FALSE)

save(arr_summary, mat_diagnostics_overall,
        file = paste0(path_aggregated_output_folder, "/aggregated_output_", 
                      model_type, ".Rdata"))

