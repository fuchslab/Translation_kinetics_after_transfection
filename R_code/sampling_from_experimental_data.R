# This R script is used to sample 100 trajectories from each of the experimental
# datasets `20160427_mean_eGFP.xlsx` and `20160427_mean_d2eGFP.xlsx`

library(readxl)

# number of observed trajectories that we want to sample from the experimental data
N <- 100

for(dataset in c("eGFP", "d2eGFP")){
  path_to_dataset <- paste0("data/experimental_data/", dataset,
                            "/20160427_mean_", dataset, ".xlsx", sep = "")
  data_mat <- as.matrix(read_excel(path_to_dataset, col_names = FALSE))
  # the first column contains the observation time points in seconds;
  # further columns contain the observed trajectories of the mean fluorescence
  # intensity of different cells

  # convert time points from hours to seconds
  time_points <- data_mat[,1]/3600


  # sample N indices of the columns containing the observed trajectories
  set.seed(2029938058)
  sampled_indices <- sort(sample.int(n = dim(data_mat)[2] - 1, size = N)) + 1

  # sample N random integers that can later be used to initialize the random
  # number generator before the MCMC sampling
  seeds <- sample.int(.Machine$integer.max, size = N)

  output_path <- paste0("data/experimental_data/", dataset, "/", sep = "")

  for(i in 1:N){
    index <- sampled_indices[i]
    random_seed <- seeds[i]
    Observations <- cbind(time = time_points, FI = data_mat[,sampled_indices[i]])
    save(Observations, index, random_seed,
         file = paste0(output_path, "observed_trajectory_", i, ".Rdata", sep = ""))
  }

  Observations <- cbind(time = time_points, FI = data_mat[,sampled_indices])
  save(Observations,
       file = paste0(output_path, "observed_trajectory_all.Rdata", sep = ""))


  # save a plot of the sampled trajectories -------------------------------------
  pdf(paste0(output_path, "sampled_trajectories.pdf", sep = ""),
      height = 6, width = 6)

  par(mar=c(3, 4, 2, 1))
  plot_color <- c('#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443',
                  '#006837','#004529')
  plot(time_points, data_mat[,sampled_indices[1]], col = plot_color[1],
       ylim = c(0, max(data_mat[,sampled_indices])),
       xlab = "", ylab = "fluorescence intensity", type = 'l', lty = 3, lwd = 3,
       las = 1)
  title(paste0("Sample of ", N, " observed trajectories", sep =""))
  for(l in 2:N){
    lines(time_points, data_mat[,sampled_indices[l]],
          col = plot_color[l %% length(plot_color) + 1], type = 'l')
  }
  mtext('time in hours', side=1, las=1, line=2)

  # in log scale
  plot(time_points, data_mat[,sampled_indices[1]], col = plot_color[1],
       ylim = c(1, max(data_mat[,sampled_indices])),
       xlab = "", ylab = "fluorescence intensity", type = 'l', lty = 3, lwd = 3,
       las = 1, log = "y")
  title(paste0("Sample of ", N, " observed trajectories", sep =""))
  for(l in 2:N){
    lines(time_points, data_mat[,sampled_indices[l]],
          col = plot_color[l %% length(plot_color) + 1], type = 'l')
  }
  mtext('time in hours', side=1, las=1, line=2)

  dev.off()
}

