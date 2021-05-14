# This Rscript is used to summarize the computational time for the sampling
# for all datasets
datasets <- c("simulated_data_dataset1_no_error", 
              "simulated_data_dataset1_with_error", 
              "experimental_data_eGFP", "experimental_data_d2eGFP")
col_names <- c("max_total_time", "mean_total_time" ) 

res_matrix_ODE <- array(rep(NA, length(datasets) * length(col_names) * 100),
                        dim= c(length(datasets), length(col_names), 100),
                        dimnames = list(datasets, col_names, 1:100))
colnames(res_matrix_ODE) <- col_names
rownames(res_matrix_ODE) <- datasets
res_matrix_SDE <- res_matrix_ODE


for( i in 1:length(datasets)){
  load(paste0("intermediate_output_files/aggregated_output/",
              datasets[i], "/aggregated_output_ODE.Rdata"))
  res_matrix_ODE[datasets[i], col_names, ] <- t(mat_diagnostics_overall[ , col_names])
  
  
  load(paste0("intermediate_output_files/aggregated_output/",
              datasets[i], "/aggregated_output_SDE.Rdata"))
  res_matrix_SDE[datasets[i], col_names, ] <- t(mat_diagnostics_overall[ , col_names])
}


mean_of_total_time <- rbind(ODE = apply(res_matrix_ODE[,"max_total_time",], 1, mean), 
                            SDE = apply(res_matrix_SDE[,"max_total_time",], 1, mean))

print("mean (over diffferent trajectories) of total sampling time:")
print(mean_of_total_time)
