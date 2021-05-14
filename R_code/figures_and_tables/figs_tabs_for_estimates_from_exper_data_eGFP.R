library(rstan)
library(RColorBrewer)

source("R_code/figures_and_tables/function_pairs_of_2stanfits.R")
source("R_code/figures_and_tables/function_Stan_output_summary.R")
source("R_code/figures_and_tables/functions_output_matrix_aggregated_results.R")
source('R_code/figures_and_tables/function_compare_traceplots.R')

save_results <- TRUE # TRUE FALSE

fig_format <- "png" # png pdf

dataset <- "experimental_data_eGFP"
obs_index <- 6 #18

# load data and stanfit objects ------------------------------------------------
# load simulated trajectory 
# split input argument `dataset` at second underscore
dataset_split <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2',
                              dataset), ' ')[[1]]


# load ODE stanfit object 
path_ODE_stanfit_object <- paste0("intermediate_output_files/stanfit_objects/",
                                  dataset , "_ODE/stanfit_object_",
                                  obs_index, ".rds", sep = "")
ODE_stanfit_object <- readRDS(file = path_ODE_stanfit_object)

# load SDE stanfit object
path_SDE_stanfit_object <- paste0("intermediate_output_files/stanfit_objects/",
                                  dataset , "_SDE/stanfit_object_",
                                  obs_index, ".rds", sep = "")
SDE_stanfit_object <- readRDS(file = path_SDE_stanfit_object)

# print output summary ODE model -----------------------------------------------
digits <- c(1,2,3,2,2,2,0,2)
param <- c("theta", "m0", "scale", "offset", "t0", "sigma", "prod_theta2_m0", 
           "prod_theta2_scale", "prod_m0_scale", "prod_theta2_m0_scale")
param_names <-c("$\\theta_1$","$\\theta_2$","$\\theta_3$", "$m_0$", 
                "$\\mathrm{scale}$", "$\\mathrm{offset}$", "$t_0$", "$\\sigma$",
                "$\\theta_2m_0$", "$\\theta_2\\mathrm{scale}$", 
                "$m_0\\mathrm{scale}$", "$\\theta_2m_0\\mathrm{scale}$")

stan_summary_ODE <- tab_Stan_output_summary(ODE_stanfit_object, param = param, 
                                            param_names = param_names)

if(save_results){
  ODE_table <- print(xtable::xtable(stan_summary_ODE, digits = digits,
                                    align = c("l", rep("r", ncol(stan_summary_ODE)))),
                     sanitize.text.function=function(x){x},
                     booktabs = TRUE,
                     file = "figures_and_tables/tab_exper_data_eGFP_ODE_summary.txt")
  ODE_table <- readLines("figures_and_tables/tab_exper_data_eGFP_ODE_summary.txt")
  writeLines(ODE_table[5:(length(ODE_table)-1)], 
             con = "figures_and_tables/tab_exper_data_eGFP_ODE_summary.txt" )
  
}else{
  print(xtable::xtable(stan_summary_ODE, digits = digits,
                       align = c("l", rep("r", ncol(stan_summary_ODE)))),
        sanitize.text.function=function(x){x},
        booktabs = TRUE)
}

# print output summary SDE model -----------------------------------------------
param <- c("theta", "m0", "scale", "offset",  "sigma", "prod_theta2_m0",
           "prod_theta2_scale", "prod_m0_scale", "prod_theta2_m0_scale")
param_names <- c("$\\theta_1$", "$\\theta_2$", "$\\theta_3$", "$m_0$", 
                 "$\\mathrm{scale}$", "$\\mathrm{offset}$", "$\\sigma$", 
                 "$\\theta_2m_0$", "$\\theta_2\\mathrm{scale}$", 
                 "$m_0\\mathrm{scale}$", "$\\theta_2m_0\\mathrm{scale}$")

stan_summary_SDE <- tab_Stan_output_summary(SDE_stanfit_object, param = param, 
                                            param_names = param_names)

if(save_results){
  SDE_table <- print(xtable::xtable(stan_summary_SDE, digits = digits,
                                    align = c("l", rep("r", ncol(stan_summary_SDE)))),
                     sanitize.text.function=function(x){x},
                     booktabs = TRUE,
                     file = "figures_and_tables/tab_exper_data_eGFP_SDE_summary.txt")
  SDE_table <- readLines("figures_and_tables/tab_exper_data_eGFP_SDE_summary.txt")
  writeLines(SDE_table[5:(length(SDE_table)-1)], 
             con = "figures_and_tables/tab_exper_data_eGFP_SDE_summary.txt" )
}else{
  print(xtable::xtable(stan_summary_SDE, digits = digits,
                       align = c("l", rep("r", ncol(stan_summary_SDE)))),
        sanitize.text.function=function(x){x},
        booktabs = TRUE)
}



# pairplots --------------------------------------------------------------------
v_col1 <- c(Helmholtz_blue, brewer.pal(9, "Blues")[9])
v_col2 <- c(Helmholtz_green, brewer.pal(9, "Greens")[c(7,9)])

if(save_results){
  if(fig_format == "pdf"){
    pdf("figures_and_tables/Fig_pairs_exper_data_eGFP_theta_1_3_offset_sigma.pdf", 
        width = 6, height = 5.3)
  }else if(fig_format == "png"){
    png("figures_and_tables/Fig_pairs_exper_data_eGFP_theta_1_3_offset_sigma.png", 
        width = 900, height = 790, pointsize = 24)
  }
}
  param <- c("theta[1]", "theta[3]", "offset", "sigma")
  pars_names <- c(expression(theta[1]), expression(theta[3]), "offset", "sigma")
  
  pairs_of_2stanfits(SDE_stanfit_object, ODE_stanfit_object,
                     pars = param, 
                     pars_names = pars_names, inc_boxplot = TRUE,
                     v_col1 = v_col1, v_col2 = v_col2,
                     name1 = expression(paste("Posterior sample for ", 
                                              bold("SDE"), " model")),
                     name2 = expression(paste("Posterior sample for ", 
                                              bold("ODE"), " model")))
if(save_results) dev.off()

if(save_results){
  if(fig_format == "pdf"){
    pdf("figures_and_tables/Fig_pairs_exper_data_eGFP_theta_2_m0_scale.pdf", 
        width = 6.35, height = 5.6)
  }else if(fig_format == "png"){
    png("figures_and_tables/Fig_pairs_exper_data_eGFP_theta_2_m0_scale.png", 
        width = 1260, height = 1140, pointsize = 32)
  }
}
  param <- c("theta[2]", "m0", "scale",
             "prod_theta2_m0", 
             "prod_theta2_m0_scale")
  pars_names <- c(expression(theta[2]), expression(m[0]), "scale", 
                  expression(theta[2]*m[0]), expression(theta[2]*m[0]*scale))
  
  pairs_of_2stanfits(SDE_stanfit_object, ODE_stanfit_object,
                     pars = param, pars_names = pars_names, inc_boxplot = TRUE,
                     v_col1 = v_col1, v_col2 = v_col2,
                     name1 = expression(paste("Posterior sample for ", 
                                              bold("SDE"), " model")),
                     name2 = expression(paste("Posterior sample for ", 
                                              bold("ODE"), " model")))
if(save_results) dev.off()

if(save_results){
  if(fig_format == "pdf"){
    pdf("figures_and_tables/Fig_pairs_exper_data_eGFP_theta_2_m0_scale_all.pdf", 
        width = 10, height = 9.4)
  }else if(fig_format == "png"){
    png("figures_and_tables/Fig_pairs_exper_data_eGFP_theta_2_m0_scale_all.png", 
        width = 2005, height = 1880, pointsize = 32)
  }
}
  param <- c("theta[2]", "m0", "scale",
             "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
             "prod_theta2_m0_scale")
  pars_names <- c(expression(theta[2]), expression(m[0]), "scale", 
                  expression(theta[2]*m[0]), expression(theta[2]*scale),
                  expression(m[0]*scale), expression(theta[2]*m[0]*scale))
  
  pairs_of_2stanfits(SDE_stanfit_object, ODE_stanfit_object,
                     pars = param, pars_names = pars_names, inc_boxplot = TRUE,
                     v_col1 = v_col1, v_col2 = v_col2,
                     name1 = expression(paste("Posterior sample for ", 
                                              bold("SDE"), " model")),
                     name2 = expression(paste("Posterior sample for ", 
                                              bold("ODE"), " model")))
if(save_results) dev.off()
  
if(save_results){
  # pdf("figures_and_tables/Fig_pairs_exper_data_eGFP_all_param.pdf", 
  #     width = 15, height = 14.1)
  png("figures_and_tables/Fig_pairs_exper_data_eGFP_all_param.png", 
      width = 2400, height = 2300, pointsize = 28)
}
param <- c("theta[1]", "theta[3]", "offset", "sigma",
           "theta[2]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale", "lp__")
pars_names <- c(expression(theta[1]), expression(theta[3]), "offset", "sigma",
                expression(theta[2]), expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*scale),
                expression(m[0]*scale), expression(theta[2]*m[0]*scale), "log-posterior")

pairs_of_2stanfits(SDE_stanfit_object, ODE_stanfit_object,
                   pars = param, pars_names = pars_names,
                   inc_boxplot = TRUE,
                   v_col1 = v_col1, v_col2 = v_col2,
                   name1 = expression(paste("Posterior sample for ", 
                                            bold("SDE"), " model")),
                   name2 = expression(paste("Posterior sample for ", 
                                            bold("ODE"), " model")),
                   incl_priors = TRUE)
if(save_results) dev.off()
  
if(save_results){
  # pdf("figures_and_tables/Fig_pairs_exper_data_eGFP_all_param_SDE.pdf", 
  #       width = 15, height = 14.1)
  png("figures_and_tables/Fig_pairs_exper_data_eGFP_all_param_SDE.png", 
      width = 2400, height = 2300, pointsize = 28)
}
param <- c("theta[1]", "theta[3]", "offset", "sigma",
           "theta[2]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale", "lp__")
pars_names <- c(expression(theta[1]), expression(theta[3]), "offset", "sigma",
                expression(theta[2]), expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*scale),
                expression(m[0]*scale), expression(theta[2]*m[0]*scale), "log-posterior")

pairs_of_2stanfits(SDE_stanfit_object,
                   pars = param, pars_names = pars_names, inc_boxplot = TRUE,
                   v_col1 = v_col1, v_col2 = v_col2,
                   name1 = expression(paste("Posterior sample for ", 
                                            bold("SDE"), " model")),
                   incl_priors = TRUE)
if(save_results) dev.off()

if(save_results){
  # pdf("figures_and_tables/Fig_pairs_exper_data_eGFP_all_param_ODE.pdf", 
  #     width = 15, height = 14.1)
  png("figures_and_tables/Fig_pairs_exper_data_eGFP_all_param_ODE.png", 
      width = 2400, height = 2300, pointsize = 28)
}
param <- c("theta[1]", "theta[3]", "offset", "sigma", 
           "theta[2]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale", "lp__")
pars_names <- c(expression(theta[1]), expression(theta[3]), "offset", "sigma",
                expression(theta[2]), expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*scale),
                expression(m[0]*scale), expression(theta[2]*m[0]*scale), "log-posterior")

pairs_of_2stanfits(ODE_stanfit_object,
                   pars = param, pars_names = pars_names, inc_boxplot = TRUE,
                   name1 = expression(paste("Posterior sample for ", 
                                            bold("ODE"), " model")),
                   col1 = Helmholtz_green,
                   v_col1 = c(Helmholtz_green, brewer.pal(9, "Blues")[9]),
                   incl_priors = TRUE)
if(save_results) dev.off()

# tables of aggregated output of all trajectories ------------------------------
# load aggregated output
load(paste0("intermediate_output_files/aggregated_output/",
            dataset, "/aggregated_output_ODE.Rdata"))
arr_summary_ODE <- arr_summary
load(paste0("intermediate_output_files/aggregated_output/",
            dataset, "/aggregated_output_SDE.Rdata"))
arr_summary_SDE <- arr_summary

param <- c("theta[1]", "theta[2]", "theta[3]", "m0", "scale",
           "prod_theta2_m0",  "prod_theta2_scale",  "prod_m0_scale",
           "prod_theta2_m0_scale", "offset", "sigma")
param_names <- c("$\\theta_1$","$\\theta_2$","$\\theta_3$", "$m_0$", 
                 "$\\mathrm{scale}$",  "$\\theta_2m_0$", 
                 "$\\theta_2\\mathrm{scale}$", "$m_0\\mathrm{scale}$",
                 "$\\theta_2m_0\\mathrm{scale}$", "$\\mathrm{offset}$",
                 "$\\sigma$")

res_all <- construct_mat_comparing_ODE_SDE(arr_summary_ODE =arr_summary_ODE, 
                                           arr_summary_SDE = arr_summary_SDE,
                                           param = param, 
                                           param_names = param_names)

if(save_results){
  file_name <- "figures_and_tables/tab_exper_data_eGFP_CI_lengths.txt"
  save_latex_table(result_table = res_all, file_name = file_name)
  
}else{
  print(xtable::xtable(res_all, digits = c(0,2,2,3), format = "latex"), 
        sanitize.text.function=function(x){x}, booktabs = TRUE)
}

# plot the observed  trajectories ----------------------------------------------
load("data/experimental_data/eGFP/observed_trajectory_all.Rdata")
Observations_eGFP <- Observations
load("data/experimental_data/d2eGFP/observed_trajectory_all.Rdata")
Observations_d2eGFP <- Observations

if(save_results){
  pdf("figures_and_tables/Fig_mRNA_transfection_data_used_in_simultion_study.pdf", height = 3.2, width = 8)
}
layout(matrix(c(1,2), ncol = 2))
par(mar=c(3.2, 4.5, 1.5, 1))
plot_color <- Helmholtz_green #c('#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')
ymax <- max(Observations_eGFP[,-1], Observations_d2eGFP[,-1])
ymin <- min(Observations_eGFP[,-1], Observations_d2eGFP[,-1])


# ---- eGFP 
plot(Observations_eGFP[,1], Observations_eGFP[ ,obs_index + 1],  
     ylim = c(ymin, ymax), las = 1,
     type = 'p', pch = 20, cex = .3, col = plot_color[1], 
     xlab = "", ylab = "")

title(main = "eGFP")
title(xlab = "time [h]", line = 2.1)
title(ylab = "fluorescence intensity [a.u.]", line = 3.5)

obs_range <- apply(Observations_eGFP[,-1], 1, range)

polygon(c(Observations_eGFP[,1], rev(Observations_eGFP[,1])), 
        c(obs_range[1,], rev(obs_range[2,])), col = "gray", density = 50)
lines(Observations_eGFP[,1], Observations_eGFP[ ,obs_index + 1],  
     type = 'p', pch = 20, cex = .3, col = plot_color[1])

# ---- d2eGFP 
plot(Observations_d2eGFP[,1], Observations_d2eGFP[ ,obs_index + 1],  ylim = c(ymin, ymax), las = 1,
     type = 'p', pch = 20, cex = .3, col = plot_color[1], 
     xlab = "", ylab = "")

title(main = "d2eGFP")
title(xlab = "time [h]", line = 2.1)
title(ylab = "fluorescence intensity [a.u.]", line = 3.5)

obs_range <- apply(Observations_d2eGFP[,-1], 1, range)

polygon(c(Observations_d2eGFP[,1], rev(Observations_d2eGFP[,1])), 
        c(obs_range[1,], rev(obs_range[2,])), col = "gray", density = 50)

if(save_results) dev.off()

# traceplots -------------------------------------------------------------------
param <- c("theta[1]", "theta[2]", "theta[3]", "m0", "scale", 
           "prod_theta2_m0", "prod_theta2_m0_scale", "offset", "sigma")
pars_names <- c(expression(theta[1]), expression(theta[2]), expression(theta[3]), 
                expression(m[0]), "scale", 
                expression(theta[2]*m[0]), expression(theta[2]*m[0]*scale),
                "offset", expression(sigma))

# create ggplot2
p <- compare_traceplots(SDE_stanfit_object, ODE_stanfit_object,
                        pars = param,
                        pars_names = pars_names,
                        inc_warmup = FALSE)


if(save_results){
  ggsave("figures_and_tables/Fig_traceplots_exper_data_eGFP.pdf", plot = p,
         width = 8, height = 9)
}else{
  p
}
