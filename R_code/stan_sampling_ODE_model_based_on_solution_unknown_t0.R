# This R script is used to perform HMC+NUTS sampling using the package rstan
# for the ODE assuming multiplicative measurement error and that the initial time
# point t0 of mRNA realease needs to be estimated.
# (For the ODE model, multiplicative measurement error is also assumed when 
# sampling for simulated data without error.)
library(rstan)

input_args <- commandArgs(trailingOnly = TRUE)

# check input arguments --------------------------------------------------------------
if(input_args[1] %in% c("experimental_data_eGFP", "experimental_data_d2eGFP",
                        "simulated_data_dataset1_with_error",
                        "simulated_data_dataset2_with_error",
                        "simulated_data_dataset1_no_error",
                        "simulated_data_dataset2_no_error")){
  dataset <- input_args[1]
} else {
  stop(paste0("Unknown dataset '",  input_args[1],
              "' supplied. Eligible first arguments are [experimental_data_eGFP,
              experimental_data_d2eGFP, simulated_data_dataset1_with_error,
              simulated_data_dataset2_with_error, simulated_data_dataset1_no_error,
              simulated_data_dataset2_no_error]."),
       sep = "")
}

obs_index <- as.numeric(input_args[2])
if(!(obs_index %in% 1:100)){
  stop("The second argument must be an integer between 1 and 100.")
}

# # load data ------------------------------------------------------------------------
# split input argument `dataset` at second underscore
dataset_split <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2',
                              dataset), ' ')[[1]]
input_path <- paste0("data/", dataset_split[1], "/", dataset_split[2], "/", sep = "")
if(dataset_split[1] == "experimental_data"){
  file_name <- paste0("observed_trajectory_", obs_index, ".Rdata", sep = "")
}else{
  file_name <- paste0("simulated_trajectory_", obs_index, ".Rdata", sep = "")
}

# contains the following variables: index, Observations, random_seed
load(file = paste0(input_path, file_name, sep = ""))

# HMC + NUTS sampling using rstan ----------------------------------------------------
rstan_options(auto_write = TRUE)
n_chains <- 8
options(mc.cores = n_chains)
n_Iter <- 5000
stan_model_file <- "Stan_model_code/ODE_model_based_on_solution_mult_error_unknown_t0.stan"

path_sample_file <- paste0("intermediate_output_files/stan_sample_files/",
                           dataset , "_ODE")
path_diagnostic_file <- paste0("intermediate_output_files/stan_diagnostic_files/",
                                dataset , "_ODE")
path_stanfit_object <- paste0("intermediate_output_files/stanfit_objects/",
                               dataset , "_ODE")
dir.create(file.path(getwd(), path_sample_file), showWarnings = FALSE)
dir.create(file.path(getwd(), path_diagnostic_file), showWarnings = FALSE)
dir.create(file.path(getwd(), path_stanfit_object), showWarnings = FALSE)

sdata <- list(
  M = dim(Observations)[1],
  y_obs = Observations[, "FI"],
  time_points = Observations[, "time"]
)

stanfit_object <-
  stan(file = stan_model_file,
       data = sdata,
       seed = random_seed,
       control = list(adapt_delta = 0.96, max_treedepth = 15),
       chains = n_chains,
       iter = n_Iter,
       verbose = TRUE,
       sample_file = paste0(path_sample_file, "/sample_file_", obs_index),
       diagnostic_file = paste0(path_diagnostic_file, "/diagnostic_file_", obs_index))


saveRDS(stanfit_object,
     file = paste0(path_stanfit_object, "/stanfit_object_", obs_index, ".rds"))

