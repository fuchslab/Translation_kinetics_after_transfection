# This R script contains several functions to construct and save a table of the 
# estimation results (in particular credible intervals) results aggregated over
# several data trajectories.

##  Define default colors ------------------------------------------------------------
HMGU_red <- '#E50030'
TUM_blue <- '#3070B3'
Helmholtz_green <- '#8CB423'
Helmholtz_blue <- '#005aa0'
Helmholtz_dark_blue <- '#0a2D6e'

## Define functions ------------------------------------------------------------
#' Construct the matrix of results aggregated over several data trajectories for 
#' one model type 
#'
#' @param arr_summary a previously aggregated 3-dim. output array
#' @param param the parameters to be included in the result matrix
#' @param true_values if available, the true parameters, default: NA
#'
#' @return  matrix of results 
construct_mat_per_stanfit <- function(arr_summary, param, true_values = NA){
  if(all(is.na(true_values))){
    res_col_names <- c("length_of_prior_hpdi", "median_length_CI", 
                       "cv_length_CI")
  }else{
    res_col_names <- c("length_of_prior_hpdi", "median_length_CI", 
                       "cv_length_CI",
                       "median_length_CI_normalized_by_true_val",
                       "num_true_value_covered")
  }
  res <- matrix(rep(NA, times = length(param) * length(res_col_names)),
                ncol = length(res_col_names))
  colnames(res) <- res_col_names
  rownames(res) <- param
  
  prior_statistics <- readRDS("intermediate_output_files/prior_statistics.rds")
  res[ , "length_of_prior_hpdi"] <- prior_statistics[param, "length of 95% HPDI"]
  
  length_CI <- arr_summary[ , , "97.5%"] - arr_summary[ , , "2.5%"]
  res[ , "median_length_CI"] <- apply(length_CI, 2, median)
  res[ , "cv_length_CI"] <- apply(length_CI, 2, function(x) var(x)/mean(x))
  # res[ , "median_length_CI_normalized_by_median"] <- 
  #   apply(length_CI / arr_summary[ , , "50%"], 2, median)
  
  if(!all(is.na(true_values))){
    true_value_in_CI <- 
      apply(arr_summary, 1,
            function(x) true_values >= x[ , "2.5%"] & true_values <= x[ , "97.5%"])
    res[ , "num_true_value_covered"] <- apply(true_value_in_CI, 1, sum)
    res[ , "median_length_CI_normalized_by_true_val"] <- 
      res[ , "median_length_CI"] / true_values
    
    # ind_tv_covered_at_least_once <- which(!is.na(res[ , "num_true_value_covered"]) & 
    #                                         res[ , "num_true_value_covered"] > 0)
    # for(i in ind_tv_covered_at_least_once){
    #   length_of_CI_covering_tv <- length_CI[true_value_in_CI[i,], i] 
    #   res[i, "median_length_CI_true_value_covered"] <-
    #     median(length_of_CI_covering_tv)
    #   res[i, "cv_length_CI_true_value_covered"] <- 
    #     var(length_of_CI_covering_tv) / mean(length_of_CI_covering_tv)
    # }
  }
  res
}


#' Construct the matrix of results aggregated over several data trajectories for 
#' the ODE and the SDE model
#'
#' @param arr_summary_ODE a previously aggregated 3-dim. output array for the ODE model
#' @param arr_summary_SDE a previously aggregated 3-dim. output array for the SDE model
#' @param param the parameters to be included in the result matrix
#' @param param_names the parameter names as to be included as row names of the 
#' summary matrix for printing, default: param
#' @param true_values if available, the true parameters, default: NA
#'
#' @return  matrix of results for both model types
construct_mat_comparing_ODE_SDE <- function(arr_summary_ODE, arr_summary_SDE,
                                            param, param_names = param, 
                                            true_values = NA){
  # select the specified parameters
  arr_summary_ODE <- arr_summary_ODE[ , param, ]  
  arr_summary_SDE <- arr_summary_SDE[ , param, ]
  # construct the results matrices for each model type separately
  res_ODE <- construct_mat_per_stanfit(arr_summary = arr_summary_ODE, 
                                       param = param, true_values = true_values)
  res_SDE <- construct_mat_per_stanfit(arr_summary = arr_summary_SDE, 
                                       param = param, true_values = true_values)
  
  # median_ratios <- apply((arr_summary_SDE[,,"97.5%"] - arr_summary_SDE[,,"2.5%"]) / 
  #                          (arr_summary_ODE[,,"97.5%"] - arr_summary_ODE[,,"2.5%"]), 2, median)
  # 
  # merge the result matrices such that every other row belongs to one model 
  # type and thus the results of the same parameter are grouped
  res_all <- rbind(res_ODE, res_SDE)
  num_row <- nrow(res_ODE)
  reorder_vec <- rep(1:num_row, each=2)
  reorder_vec[2*(1:num_row)] <- reorder_vec[2*(1:num_row)] + num_row
  res_all <- res_all[reorder_vec,]
  rownames(res_all) <- 
    paste(rep(param_names, each =2), 
          rep(c("ODE", "SDE"), times = nrow(res_ODE)), sep =" & ")
  res_all
}


#' Saves the results table to .txt-file in latex format
#'
#' @param result_table matrix of results
#' @param digits vector of the number of digits to be printed per column; 
#' there is an extra 0 at the beginning for the column of the row names
#' @param file_name name of the file to which the table is to be saved
#' @param inc_true_values boolean indicating whether results depending on the 
#' true parameter values are included
#'
#' @return 
save_latex_table <- function(result_table, digits = NA, file_name,
                             inc_true_values = FALSE){
  if(all(is.na(digits))){
    if(inc_true_values){
      digits <- c(0,2,2,3,2,0)
    }else{
      digits <- c(0,2,2,3)
    }
  }
  print(xtable::xtable(result_table, digits = digits, format = "latex"), 
        sanitize.text.function=function(x){x}, booktabs = TRUE,
        file = file_name)

  res_table <- readLines(file_name)
  res_table_content <- res_table[9:(length(res_table)-3)]
  for(i in 1:length(res_table_content)){
    split <- stringr::str_split(res_table_content[i], '&', n=2)[[1]]
    if(i %% 2 == 1){
      res_table_content[i] <- paste0("  \\multirow{2}{*}{", split[1], "} & ", split[2])
    }else{
      if(i < length(res_table_content)){
        res_table_content[i] <- paste0("    & ", split[2], "\\hline")
      }else{
        res_table_content[i] <- paste0("    & ", split[2])
      }
      
    }
  }
  if(inc_true_values){
    header <- 
      c("\\begin{tabular}{llrrrrrr}",
        "  \\toprule",
        "  & & length of      & median      & c.v. of  		& median of     & number \\\\ ",
        "  & & prior  95\\%   & length of   & lengths of  & length of CIs & of CIs \\\\ ", 
        "  & & center         & 95\\% CIs   & 95\\% CIs   & rescaled by   & covering\\\\ ",
        "  & & interval       &             &             & true value    & true value\\\\ ",
        "  \\midrule")
      # c("\\begin{tabular}{llrrrrr}",
      #   "  \\toprule",
      #   "  & & median         & c.v. of  		  & number of     & median length   &  c.v. of  lengths\\\\ ",
      #   "  & & length         & lengths 		  & CIs covering  & of CIs covering &  of CIs covering \\\\ ", 
      #   "  & & of 95\\% CIs   & of 95\\% CIs  & true value    & true value      &  true value \\\\ ",
      #   "  \\midrule")
  }else{
    header <- 
      c("\\begin{tabular}{llrrrrr}",
        "  \\toprule",
        "  & & length of    & median      & c.v. of  		\\\\ ",
        "  & & prior 95\\%  & length of   & lengths of	\\\\ ", 
        "  & & center       & 95\\% CIs   & 95\\% CIs   \\\\ ",
        "  & & interval     &             &             \\\\ ",
        "  \\midrule")
  }

  bottom <- c("  \\bottomrule", "\\end{tabular}")
  res_table <- c(header, res_table_content, bottom)
  writeLines(res_table, con = file_name )
  invisible(res_table)
}


plot_aggregate_values <- function(arr_summary_1, arr_summary_2,
                                  name1 = "ODE", name2 = "SDE",
                                  param, 
                                  param_names = param,
                                  y_values_name = "median_length_CI_normalized_by_true_val",
                                  ylab = "median of lengths of CIs rescaled by true value",
                                  y_log_scale = FALSE, 
                                  true_values = NA){
  layout(matrix(1, ncol = 1))
  # select the specified parameters
  arr_summary_1 <- arr_summary_1[ , param, ]  
  arr_summary_2 <- arr_summary_2[ , param, ]
  # construct the results matrices for each model type separately
  res_ODE <- construct_mat_per_stanfit(arr_summary = arr_summary_1, 
                                       param = param, true_values = true_values)
  res_SDE <- construct_mat_per_stanfit(arr_summary = arr_summary_2, 
                                       param = param, true_values = true_values)
  
  if(y_log_scale){
    log_scale <- "y"
  }else{
    log_scale <- ""
  }
  
  col <- c('#67001f', HMGU_red,'#d6604d','#f4a582','#fddbc7','#d1e5f0',
           '#92c5de','#4393c3', TUM_blue, Helmholtz_dark_blue, '#081d58', 1)
  ylim <- range(res_ODE[ ,y_values_name], res_SDE[ ,y_values_name])
  
  par(xpd=TRUE, mar = c(3.2,4.5,2,7))
  plot(res_ODE[ ,"num_true_value_covered"], 
       res_ODE[ ,y_values_name],
       pch = 1, las = 1, col=col,
       log = log_scale,
       xlab = "",
       ylab = "",
       xlim = c(0,100),
       ylim = ylim)
  title(xlab = "number of CIs covering the true value", line = 2.2)
  title(ylab = ylab, line = 3.5)
  lines(res_SDE[,"num_true_value_covered"], 
        res_SDE[,y_values_name],
        type = 'p', pch = 2, col = col)
  
  legend("topright", inset=c(-.385,0), legend = c(name1, name2, param_names),
         pch = c(1, 2, rep(15, times = length(param))), 
         col = c(1,1 , col[1:length(param)]), bty='n')
  
}