library(rstan)
library(RColorBrewer)

source("R_code/figures_and_tables/function_pairs_of_2stanfits.R")
source("R_code/figures_and_tables/functions_output_matrix_aggregated_results.R")

save_results <- FALSE # TRUE FALSE

fig_format <- "png" # png pdf

dataset <- "simulated_data_dataset1"
obs_index <- 6

# load data and stanfit objects ------------------------------------------------
# load simulated trajectory 
# split input argument `dataset` at second underscore
dataset_split <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2',
                              dataset), ' ')[[1]]
load(paste0("data/", dataset_split[1], "/", dataset_split[2], 
            "_with_error/simulated_trajectory_", obs_index, ".Rdata"))

# load ODE stanfit object for simulated data with error
path_ODE_stanfit_object_with_error <- paste0("intermediate_output_files/stanfit_objects/",
                                  dataset , "_with_error_ODE/stanfit_object_",
                                  obs_index, ".rds", sep = "")
ODE_stanfit_object_with_error <- readRDS(file = path_ODE_stanfit_object_with_error)

# load SDE stanfit object for simulated data with error
path_SDE_stanfit_object_with_error <- paste0("intermediate_output_files/stanfit_objects/",
                                  dataset , "_with_error_SDE/stanfit_object_",
                                  obs_index, ".rds", sep = "")
SDE_stanfit_object_with_error <- readRDS(file = path_SDE_stanfit_object_with_error)

# load ODE stanfit object for simulated data without error
path_ODE_stanfit_object_no_error <- paste0("intermediate_output_files/stanfit_objects/",
                                             dataset , "_no_error_ODE/stanfit_object_",
                                             obs_index, ".rds", sep = "")
ODE_stanfit_object_no_error <- readRDS(file = path_ODE_stanfit_object_no_error)

# load SDE stanfit object for simulated data without error
path_SDE_stanfit_object_no_error <- paste0("intermediate_output_files/stanfit_objects/",
                                             dataset , "_no_error_SDE/stanfit_object_",
                                             obs_index, ".rds", sep = "")
SDE_stanfit_object_no_error <- readRDS(file = path_SDE_stanfit_object_no_error)


# pairplots to compare ODE results for data with and without error---------------------
v_col1 <- c(Helmholtz_green, brewer.pal(9, "Blues")[9])
v_col2 <- c(HMGU_red, brewer.pal(9, "Reds")[c(9)])
col1 <- Helmholtz_green
col2 <- HMGU_red
if(save_results){
  if(fig_format == "pdf"){
    pdf("figures_and_tables/Fig_pairs_ODE_sim_data_no_and_with_error_theta_1_3_offset_t0_sigma.pdf", 
        width = 7, height = 6.3)
  }else if(fig_format == "png"){
    png("figures_and_tables/Fig_pairs_ODE_sim_data_no_and_with_error_theta_1_3_offset_t0_sigma.png", 
        width = 1060, height = 960, pointsize = 24)
  }
}
param <- c("theta[1]", "theta[3]", "offset", "t0", "sigma")
true_values <- c(pars$theta[c(1,3)], pars$offset, pars$t0, pars$sigma)
pars_names <- c(expression(theta[1]), expression(theta[3]), "offset",
                expression(t[0]), expression(sigma))

pairs_of_2stanfits(ODE_stanfit_object_no_error, ODE_stanfit_object_with_error, 
                   name1 = expression(paste("Posterior sample for ODE model given simulated data ", bold("without"), " error")),
                   name2 = expression(paste("Posterior sample for ODE model given simulated data ", bold("with"), " error")),
                   pars = param, 
                   pars_names = pars_names,
                   true_values = true_values, inc_boxplot = TRUE,
                   col1 = col1, col2 = col2,
                   v_col1 = v_col1, v_col2 = v_col2)

if(save_results) dev.off()

if(save_results){
  if(fig_format == "pdf"){
    pdf("figures_and_tables/Fig_pairs_ODE_sim_data_no_and_with_error_theta_2_m0_scale.pdf", 
        width = 6.35, height = 5.6)
  }else if(fig_format == "png"){
    png("figures_and_tables/Fig_pairs_ODE_sim_data_no_and_with_error_theta_2_m0_scale.png", 
        width = 1260, height = 1140, pointsize = 32)
  }
}
param <- c("theta[2]", "m0", "scale",
           "prod_theta2_m0", "prod_theta2_m0_scale")
true_values <- c(pars$theta[2], pars$m0, pars$scale,
                 pars$theta[2]*pars$m0, pars$theta[2]*pars$m0*pars$scale)
pars_names <- c(expression(theta[2]), expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*m[0]*scale))

pairs_of_2stanfits(ODE_stanfit_object_no_error, ODE_stanfit_object_with_error,
                   name1 = expression(paste("Posterior sample for ODE model given simulated data ", bold("without"), " error")),
                   name2 = expression(paste("Posterior sample for ODE model given simulated data ", bold("with"), " error")),
                   pars = param, pars_names = pars_names,
                   true_values = true_values, inc_boxplot = TRUE,
                   col1 = col1, col2 = col2,
                   v_col1 = v_col1, v_col2 = v_col2)

if(save_results) dev.off()


if(save_results){
  if(fig_format == "pdf"){
    pdf("figures_and_tables/Fig_pairs_ODE_sim_data_no_and_with_error_theta_2_m0_scale_all.pdf", 
        width = 8.5, height = 8)
  }else if(fig_format == "png"){
    png("figures_and_tables/Fig_pairs_ODE_sim_data_no_and_with_error_theta_2_m0_scale_all.png", 
        width = 1700, height = 1600, pointsize = 32)
  }
} 
param <- c("theta[2]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale")
true_values <- c(pars$theta[2], pars$m0, pars$scale,
                 pars$theta[2]*pars$m0, pars$theta[2]*pars$scale, 
                 pars$scale*pars$m0, pars$theta[2]*pars$m0*pars$scale)
pars_names <- c(expression(theta[2]), expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*scale),
                expression(m[0]*scale), expression(theta[2]*m[0]*scale))

pairs_of_2stanfits(ODE_stanfit_object_no_error, ODE_stanfit_object_with_error,
                   name1 = expression(paste("Posterior sample for ODE model given simulated data ", bold("without"), " error")),
                   name2 = expression(paste("Posterior sample for ODE model given simulated data ", bold("with"), " error")),
                   pars = param, pars_names = pars_names,
                   true_values = true_values, inc_boxplot = TRUE,
                   col1 = col1, col2 = col2,
                   v_col1 = v_col1, v_col2 = v_col2)

if(save_results) dev.off()


if(save_results){
  # pdf("figures_and_tables/Fig_pairs_sim_data_no_and_with_error_all_param_ODE.pdf", 
  #     width = 15, height = 14.1)
  png("figures_and_tables/Fig_pairs_sim_data_no_and_with_error_all_param_ODE.png", 
      width = 2400, height = 2300, pointsize = 28)
}
param <- c("theta[1]", "theta[3]", "offset", "t0", "sigma", 
           "theta[2]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale", "lp__")
true_values <- c(pars$theta[c(1,3)], pars$offset, pars$t0, pars$sigma,
                 pars$theta[2], pars$m0, pars$scale,
                 pars$theta[2]*pars$m0, pars$theta[2]*pars$scale, 
                 pars$scale*pars$m0, pars$theta[2]*pars$m0*pars$scale, NA)
pars_names <- c(expression(theta[1]), expression(theta[3]), "offset", 
                expression(t[0]), "sigma",
                expression(theta[2]), expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*scale),
                expression(m[0]*scale), expression(theta[2]*m[0]*scale), "log-posterior")

pairs_of_2stanfits(ODE_stanfit_object_no_error, ODE_stanfit_object_with_error,
                   name1 = expression(paste("Posterior sample for ODE model given simulated data ", bold("without"), " error")),
                   name2 = expression(paste("Posterior sample for ODE model given simulated data ", bold("with"), " error")),
                   pars = param, pars_names = pars_names,
                   true_values = true_values, inc_boxplot = TRUE,
                   col1 = col1, col2 = col2,
                   v_col1 = v_col1, v_col2 = v_col2,
                   incl_priors = TRUE)
if(save_results) dev.off()



# pairplots to compare SDE results for data with and without error--------------
v_col1 <- c(Helmholtz_blue, brewer.pal(9, "Blues")[9])
v_col2 <- c(HMGU_red, brewer.pal(9, "Reds")[c(9)])
col1 <- Helmholtz_blue
col2 <- HMGU_red
if(save_results){
  if(fig_format == "pdf"){
    pdf("figures_and_tables/Fig_pairs_SDE_sim_data_no_and_with_error_theta_1_3.pdf", 
        width = 6, height = 5.3)
  }else if(fig_format == "png"){
    png("figures_and_tables/Fig_pairs_SDE_sim_data_no_and_with_error_theta_1_3.png", 
        width = 900, height = 790, pointsize = 24)
  }
}
param <- c("theta[1]", "theta[3]")
true_values <- c(pars$theta[c(1,3)])
pars_names <- c(expression(theta[1]), expression(theta[3]))

pairs_of_2stanfits(SDE_stanfit_object_no_error, SDE_stanfit_object_with_error, 
                   name1 = expression(paste("Posterior sample for SDE model given simulated data ", bold("without"), " error")),
                   name2 = expression(paste("Posterior sample for SDE model given simulated data ", bold("with"), " error")),
                   pars = param, 
                   pars_names = pars_names,
                   true_values = true_values, inc_boxplot = TRUE,
                   col1 = col1, col2 = col2,
                   v_col1 = v_col1, v_col2 = v_col2,
                   cex_title = 1)

if(save_results) dev.off()

if(save_results){
  if(fig_format == "pdf"){
    pdf("figures_and_tables/Fig_pairs_SDE_sim_data_no_and_with_error_theta_2_m0_scale.pdf", 
        width = 6.35, height = 5.6)
  }else if(fig_format == "png"){
    png("figures_and_tables/Fig_pairs_SDE_sim_data_no_and_with_error_theta_2_m0_scale.png", 
        width = 1260, height = 1140, pointsize = 32)
  }
}
param <- c("theta[2]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_m0_scale")
true_values <- c(pars$theta[2], pars$m0, pars$scale,
                 pars$theta[2]*pars$m0, pars$theta[2]*pars$m0*pars$scale)
pars_names <- c(expression(theta[2]), expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*m[0]*scale))

pairs_of_2stanfits(SDE_stanfit_object_no_error, SDE_stanfit_object_with_error, 
                   name1 = expression(paste("Posterior sample for SDE model given simulated data ", bold("without"), " error")),
                   name2 = expression(paste("Posterior sample for SDE model given simulated data ", bold("with"), " error")),
                   pars = param, pars_names = pars_names,
                   true_values = true_values, inc_boxplot = TRUE,
                   col1 = col1, col2 = col2,
                   v_col1 = v_col1, v_col2 = v_col2)

if(save_results) dev.off()

if(save_results){
  if(fig_format == "pdf"){
    pdf("figures_and_tables/Fig_pairs_SDE_sim_data_no_and_with_error_theta_2_m0_scale_all.pdf", 
        width = 8.5, height = 8)
  }else if(fig_format == "png"){
    png("figures_and_tables/Fig_pairs_SDE_sim_data_no_and_with_error_theta_2_m0_scale_all.png", 
        width = 1700, height = 1600, pointsize = 32)
  }
}
param <- c("theta[2]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale")
true_values <- c(pars$theta[2], pars$m0, pars$scale,
                 pars$theta[2]*pars$m0, pars$theta[2]*pars$scale, 
                 pars$scale*pars$m0, pars$theta[2]*pars$m0*pars$scale)
pars_names <- c(expression(theta[2]), expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*scale),
                expression(m[0]*scale), expression(theta[2]*m[0]*scale))

pairs_of_2stanfits(SDE_stanfit_object_no_error, SDE_stanfit_object_with_error, 
                   name1 = expression(paste("Posterior sample for SDE model given simulated data ", bold("without"), " error")),
                   name2 = expression(paste("Posterior sample for SDE model given simulated data ", bold("with"), " error")),
                   pars = param, pars_names = pars_names,
                   true_values = true_values, inc_boxplot = TRUE,
                   col1 = col1, col2 = col2,
                   v_col1 = v_col1, v_col2 = v_col2)

if(save_results) dev.off()

if(save_results){
  # pdf("figures_and_tables/Fig_pairs_sim_data_no_and_with_error_all_param_SDE.pdf", 
  #     width = 15, height = 14.1)
  png("figures_and_tables/Fig_pairs_sim_data_no_and_with_error_all_param_SDE.png", 
      width = 2400, height = 2300, pointsize = 28)
}
param <- c("theta[1]", "theta[3]", 
           "theta[2]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale", "lp__")
true_values <- c(pars$theta[c(1,3)],
                 pars$theta[2], pars$m0, pars$scale,
                 pars$theta[2]*pars$m0, pars$theta[2]*pars$scale, 
                 pars$scale*pars$m0, pars$theta[2]*pars$m0*pars$scale, NA)
pars_names <- c(expression(theta[1]), expression(theta[3]), 
                expression(theta[2]), expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*scale),
                expression(m[0]*scale), expression(theta[2]*m[0]*scale), "log-posterior")

pairs_of_2stanfits(SDE_stanfit_object_no_error, SDE_stanfit_object_with_error,
                   name1 = expression(paste("Posterior sample for SDE model given simulated data ", bold("without"), " error")),
                   name2 = expression(paste("Posterior sample for SDE model given simulated data ", bold("with"), " error")),
                   pars = param, pars_names = pars_names,
                   true_values = true_values, inc_boxplot = TRUE,
                   col1 = col1, col2 = col2,
                   v_col1 = v_col1, v_col2 = v_col2,
                   incl_priors = TRUE)
if(save_results) dev.off()

# plot of aggregated sampling output ODE-----------------------------------------------
# load aggregated output
load(paste0("intermediate_output_files/aggregated_output/",
            dataset, "_no_error/aggregated_output_ODE.Rdata"))
arr_summary_ODE_no_err <- arr_summary
load(paste0("intermediate_output_files/aggregated_output/",
            dataset, "_with_error/aggregated_output_ODE.Rdata"))
arr_summary_ODE_with_err <- arr_summary
param <- c("theta[1]", "theta[2]", "theta[3]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale", "offset", "sigma", "t0")
param_names <- c(expression(theta[1]), expression(theta[2]), expression(theta[3]), 
                 expression(m[0]), "scale", 
                 expression(theta[2]*m[0]), expression(theta[2]*scale),
                 expression(m[0]*scale), expression(theta[2]*m[0]*scale),
                 "offset", expression(sigma), expression(t[0]))
true_values <- c(pars$theta[1], pars$theta[2], pars$theta[3], pars$m0, pars$scale,
                 pars$theta[2]*pars$m0, pars$theta[2]*pars$scale, 
                 pars$scale*pars$m0, pars$theta[2]*pars$m0*pars$scale,
                 pars$offset, pars$sigma, pars$t0)

if(save_results){
  pdf("figures_and_tables/Fig_aggregated_sampling_output_sim_data_ODE_model_median_of_CI_length_normalized_by_true_value.pdf", 
      width = 5.7, height = 4.4)
}
plot_aggregate_values(arr_summary_1 = arr_summary_ODE_no_err, 
                      arr_summary_2 = arr_summary_ODE_with_err,
                      name1 = "without error", name2 = "with error",
                      param = param, 
                      param_names = param_names,
                      y_values_name = "median_length_CI_normalized_by_true_val",
                      ylab = "median of lengths of CIs rescaled by true value",
                      y_log_scale = TRUE,
                      true_values = true_values)

if(save_results) dev.off()

# plot of aggregated sampling output SDE-----------------------------------------------
# load aggregated output
load(paste0("intermediate_output_files/aggregated_output/",
            dataset, "_no_error/aggregated_output_SDE.Rdata"))
arr_summary_SDE_no_err <- arr_summary
load(paste0("intermediate_output_files/aggregated_output/",
            dataset, "_with_error/aggregated_output_SDE.Rdata"))
arr_summary_SDE_with_err <- arr_summary
param <- c("theta[1]", "theta[2]", "theta[3]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale")
param_names <- c(expression(theta[1]), expression(theta[2]), expression(theta[3]), 
                 expression(m[0]), "scale", 
                 expression(theta[2]*m[0]), expression(theta[2]*scale),
                 expression(m[0]*scale), expression(theta[2]*m[0]*scale))
true_values <- c(pars$theta[1], pars$theta[2], pars$theta[3], pars$m0, pars$scale,
                 pars$theta[2]*pars$m0, pars$theta[2]*pars$scale, 
                 pars$scale*pars$m0, pars$theta[2]*pars$m0*pars$scale)

if(save_results){
  pdf("figures_and_tables/Fig_aggregated_sampling_output_sim_data_SDE_model_median_of_CI_length_normalized_by_true_value.pdf", 
      width = 5.7, height = 4.4)
}
plot_aggregate_values(arr_summary_1 = arr_summary_SDE_no_err, 
                      arr_summary_2 = arr_summary_SDE_with_err,
                      name1 = "without error", name2 = "with error",
                      param = param, 
                      param_names = param_names,
                      y_values_name = "median_length_CI_normalized_by_true_val",
                      ylab = "median of lengths of CIs rescaled by true value",
                      y_log_scale = TRUE,
                      true_values = true_values)

if(save_results) dev.off()



