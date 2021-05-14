


datasets <- c("simulated_data_dataset1_no_error", 
              "simulated_data_dataset1_with_error", 
              "experimental_data_eGFP", "experimental_data_d2eGFP")
# pre-allocate memory for summary matrix of the Stan output diagnostics
# that take integer values (num. div. transitions and num. treedepth exceeded)
col_names <- c("no", "1-10", "11-100", ">100", "maximum") 
res_matrix_ODE_div_trans <- matrix(rep(NA, length(datasets) * length(col_names)),
                                   nrow = length(datasets), ncol = length(col_names),
                                   dimnames = list(datasets, col_names))
res_matrix_SDE_div_trans <- res_matrix_SDE_max_td_exc <- 
  res_matrix_ODE_max_td_exc <- res_matrix_ODE_div_trans

# function to calculate the rows of the summary matrix of the Stan output diagnostics
# that take integer values
calc_row_int_diag <- function(mat_diagnostics_overall, colname){
  col_of_interest <- mat_diagnostics_overall[, colname]
  mat_row <- c(sum(col_of_interest == 0), 
               sum(col_of_interest > 0 & col_of_interest <= 10),
               sum(col_of_interest > 10 & col_of_interest <= 100),
               sum(col_of_interest > 100),
               max(col_of_interest))
  mat_row
}

# pre-allocate memory for summary matrix of the Stan output diagnostics BFMI
col_names_bfmi <- c("mean of minima", "s.d. of minima", "mean of means", "s.d. of means") 
res_matrix_ODE_bfmi <- matrix(rep(NA, length(datasets) * length(col_names_bfmi)),
                              nrow = length(datasets), ncol = length(col_names_bfmi),
                              dimnames = list(datasets, col_names_bfmi))
res_matrix_SDE_bfmi <- res_matrix_ODE_bfmi

# function to calculate the rows of the summary matrix of the Stan output diagnostics
# that take integer values
calc_row_bfmi <- function(mat_diagnostics_overall){
  col_min_bfmi <- mat_diagnostics_overall[, "min_BFMI"]
  col_mean_bfmi <- mat_diagnostics_overall[, "mean_BFMI"]
  mat_row <- c(mean(col_min_bfmi), sd(col_min_bfmi),
               mean(col_mean_bfmi), sd(col_mean_bfmi))
  mat_row
}


for( i in 1:length(datasets)){
  # compile summary table for ODE
  load(paste0("intermediate_output_files/aggregated_output/",
              datasets[i], "/aggregated_output_ODE.Rdata"))
  res_matrix_ODE_div_trans[datasets[i], ] <- 
    calc_row_int_diag(mat_diagnostics_overall, "num_div_trans")
  res_matrix_ODE_max_td_exc[datasets[i], ] <- 
    calc_row_int_diag(mat_diagnostics_overall, "num_treedepth_exc")
  res_matrix_ODE_bfmi[datasets[i], ] <- calc_row_bfmi(mat_diagnostics_overall)
  
  # compile summary table for SDE
  load(paste0("intermediate_output_files/aggregated_output/",
              datasets[i], "/aggregated_output_SDE.Rdata"))
  res_matrix_SDE_div_trans[datasets[i], ] <- 
    calc_row_int_diag(mat_diagnostics_overall, "num_div_trans")
  res_matrix_SDE_max_td_exc[datasets[i], ] <- 
    calc_row_int_diag(mat_diagnostics_overall, "num_treedepth_exc")
  res_matrix_SDE_bfmi[datasets[i], ] <- calc_row_bfmi(mat_diagnostics_overall)
}

# save latex tables ------------------------------------------------------------
row_names <- c("simulated data without error", 
              "simulated data with error", 
              "experimental data for eGFP", 
              "experimental data for d2eGFP")

save_latex_table <- function(result_table, digits = NA, file_name, row_names,
                             table_type){
  rownames(result_table) <- row_names
  print(xtable::xtable(result_table, digits = digits, format = "latex"), 
        sanitize.text.function=function(x){x}, booktabs = TRUE,
        file = file_name)
  
  res_table <- readLines(file_name)
  res_table_content <- res_table[9:(length(res_table)-3)]
  if(table_type == "bfmi"){
    header <- c("\\begin{tabular}{lrrrr}",
                "  \\toprule",
                "  \\multirow{2}{*}{dataset} &  mean of & s.d. of &  mean of & s.d. of \\\\ ",
                "                            &  minima  & minima  &  means	  & means  \\\\ ", 
                "  \\midrule")
  }else{
    header <- c("\\begin{tabular}{lrrrrr}",
                "  \\toprule",
                "  dataset & none &  $1 - 10$ & $11 - 100$ & $> 100$ & maximum \\\\ ",
                "  \\midrule")
  }
  bottom <- c("  \\bottomrule", "\\end{tabular}")
  res_table <- c(header, res_table_content, bottom)
  writeLines(res_table, con = file_name )
  invisible(res_table)
}


digits <- rep(0, times = ncol(res_matrix_SDE_div_trans) + 1)
save_latex_table(result_table = res_matrix_SDE_div_trans, digits, 
                 file_name = "figures_and_tables/tab_Stan_diag_divergent_transitions_SDE.txt", 
                 row_names, table_type = "div_trans")
save_latex_table(result_table = res_matrix_ODE_div_trans, digits, 
                 file_name = "figures_and_tables/tab_Stan_diag_divergent_transitions_ODE.txt", 
                 row_names, table_type = "div_trans")

save_latex_table(result_table = res_matrix_SDE_max_td_exc, digits, 
                 file_name = "figures_and_tables/tab_Stan_diag_max_treedepth_exceeded_SDE.txt", 
                 row_names, table_type = "max_treedepth_exceeded")
save_latex_table(result_table = res_matrix_ODE_max_td_exc, digits, 
                 file_name = "figures_and_tables/tab_Stan_diag_max_treedepth_exceeded_ODE.txt", 
                 row_names, table_type = "max_treedepth_exceeded")

digits <- rep(2, times = ncol(res_matrix_SDE_bfmi) + 1)
save_latex_table(result_table = res_matrix_SDE_bfmi, digits, 
                 file_name = "figures_and_tables/tab_Stan_diag_BFMI_SDE.txt", 
                 row_names, table_type = "bfmi")
save_latex_table(result_table = res_matrix_ODE_bfmi, digits, 
                 file_name = "figures_and_tables/tab_Stan_diag_BFMI_ODE.txt", 
                 row_names, table_type = "bfmi")