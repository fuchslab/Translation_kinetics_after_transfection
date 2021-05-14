# This R script is used to simulate 100 trajectories with Gillespie's algorithm
# for given parameter combinations.

source("R_code/functions_to_simulate_dynamical_systems/create_model.R")
source("R_code/functions_to_simulate_dynamical_systems/gillespie.R")

save_results <- TRUE

input_args <- commandArgs(trailingOnly = TRUE)

# # check input arguments --------------------------------------------------------------
if(input_args[1] %in% c("dataset1_no_error", "dataset2_no_error",
                        "dataset1_with_error", "dataset2_with_error")){
  dataset <- input_args[1]
} else {
  stop(paste0("Unknown dataset '",  input_args[1],
              "' supplied. Eligible first arguments are [dataset1_no_error,
              dataset1_with_error, dataset2_no_error, dataset2_with_error]."),
       sep = "")
}
#dataset <- "dataset2_with_error"

# create output path and corresponding directory
output_path <- paste0("data/simulated_data/", dataset, "/", sep = "")
dir.create(file.path(getwd(), output_path), showWarnings = FALSE)

# define the model parameters
if(dataset == "dataset1_no_error"){
  pars <- list(theta = c(0.2, 0.32, 0.01), m0 = 240, t0 = 0.96, 
               scale = 1.8, offset = 6.5, sigma = 0)
}else if(dataset == "dataset1_with_error"){
  pars <- list(theta = c(0.2, 0.32, 0.01), m0 = 240, t0 = 0.96, 
               scale = 1.8, offset = 6.5, sigma = 0.02)
}else if(dataset == "dataset2_no_error"){
  pars <- list(theta = c(0.16, 0.6, 0.09), m0 = 30, t0 = 0.58, 
               scale = 1.6, offset = 9, sigma = 0)
}else if(dataset == "dataset2_with_error"){
  pars <- list(theta = c(0.16, 0.6, 0.09), m0 = 30, t0 = 0.58, 
               scale = 1.6, offset = 9, sigma = 0.05)
} 

# number of paths to simulate
num_paths <- 100

# sample random integers that can later be used to initialize the random
# number generator before the MCMC sampling
set.seed(1373268273)
seeds <- sample.int(.Machine$integer.max, size = num_paths)

# mRNA transfection - model 1
# initial conditions for mRNA, GFP
initial <- c(pars$m0, 0)
# parameters:
kinetic_pars <- pars$theta
# now create a model object:
model <- create_model_mRNA1(initial = initial, theta = kinetic_pars)

dim <- length(model$initial)

## Simulation - shifting by initial time point and taking observations ---------
# at equidistant time steps 
T <- 30 # maximum time point
dt <- 1/6 # time step
M <- T / dt + 1 # total number of observed time points
# pre-allocate memory for all trajectories
gill_process_states <- array(dim = c(M , dim + 1, num_paths),
             dimnames = list(time = NULL, 
                             path_components = c("time", "mRNA", "GFP"),
                             path_nr = NULL))
# simulate each trajectory individually
for(i in 1:num_paths){
  # simulate one trajectory with Gillespie's algorithm (don't use equidistant 
  # option of the function because of time shift; and simulation speed)
  gill_process_states1 <- gillespie(model, maxtime = T)
  
  # time shift trajectory and take time equidistant observations
  gill_process_states1[ ,1] <- gill_process_states1[ ,1] + pars$t0
  gill_process_states1 <- rbind(vector('numeric', length = dim + 1), 
                                gill_process_states1)
  get_equidist_time_shifted_obs <- function(xmat, maxtime, dt){
    numCol <- ncol(xmat)
    time <- seq(from = 0, to = maxtime, by = dt)
    indizes <- sapply(time, function(x) max(which(xmat[,1]<=x)))
    res_mat <- matrix(NA, ncol = numCol, nrow = length(indizes))
    res_mat[, 2:numCol] <- xmat[indizes, 2:numCol]
    res_mat[, 1] <- time
    res_mat
  }
  gill_process_states[ , , i] <- 
    get_equidist_time_shifted_obs(gill_process_states1, maxtime = T, dt = dt)
}

# obtain observations by linear transformation (and adding measurement error)----
obs_scale_mult_err <- gill_process_states[ ,3, ]  * pars$scale + pars$offset
obs_scale_mult_err <- obs_scale_mult_err * 
  matrix(rlnorm( n = M * num_paths, mean = 0, sd = pars$sigma),
        ncol = num_paths)
  
 

# save individual observations -------------------------------------------------
time_points <- gill_process_states[ ,1,1]
data_mat <- obs_scale_mult_err

if(save_results){
  
  for(i in 1:num_paths){
    random_seed <- seeds[i]
    Observations <- cbind(time = time_points, FI = data_mat[,i])
    Process_states <- gill_process_states[ , , i]
    save(Observations, random_seed, pars, Process_states,
         file = paste0(output_path, "simulated_trajectory_", i, ".Rdata", sep = ""))
  }

}
  
# save a plot of the sampled trajectories --------------------------------------
if(save_results){
  pdf(paste0(output_path, "simulated_trajectories.pdf", sep = ""),
    height = 6, width = 6)
}

  par(mar=c(3, 4, 3, 1))
  plot_color <- c('#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443',
                  '#006837','#004529')
  plot(time_points, data_mat[,1], col = plot_color[1], 
       ylim = c(0, max(data_mat)),
       xlab = "", ylab = "fluorescence intensity", type = 'l',
       las = 1)
  title(paste0(num_paths, " simulated trajectories with theta = (", 
               paste0(pars$theta,collapse= ", "), "), \nm0 = ", pars$m0,
               ", t0 = ", pars$t0, ", scale = ", pars$scale, ", offset = ", 
               pars$offset, ", sigma = ", pars$sigma, sep =""))
  for(l in 2:num_paths){
    lines(time_points, data_mat[,l],
          col = plot_color[l %% length(plot_color) + 1], type = 'l')
  }
  mtext('time in hours', side=1, las=1, line=2)
  
  # in log scale
  plot(time_points, data_mat[,1], col = plot_color[1],
       ylim = c(1, max(data_mat)),
       xlab = "", ylab = "fluorescence intensity", type = 'l',
       las = 1, log = "y")
  title(paste0(num_paths, " simulated trajectories with theta = (", 
               paste0(pars$theta,collapse= ", "), "), \nm0 = ", pars$m0,
               ", t0 = ", pars$t0, ", scale = ", pars$scale, ", offset = ", 
               pars$offset, ", sigma = ", pars$sigma, sep =""))
  for(l in 2:num_paths){
    lines(time_points, data_mat[,l],
          col = plot_color[l %% length(plot_color) + 1], type = 'l')
  }
  mtext('time in hours', side=1, las=1, line=2)

if(save_results) dev.off()
