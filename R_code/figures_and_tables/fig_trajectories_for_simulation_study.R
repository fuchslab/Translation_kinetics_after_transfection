# This R script is used to generate a plot of the data used in the simulation study.  

save_plots <- TRUE # TRUE FALSE


Helmholtz_green <- '#8CB423'
Helmholtz_blue <- '#005aa0'
Helmholtz_dark_blue <- '#0a2D6e'

dataset <- "simulated_data_dataset1"
obs_index <- 6

# load simulated trajectory
# split input argument `dataset` at second underscore
dataset_split <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2',
                              dataset), ' ')[[1]]

load(paste0("data/", dataset_split[1], "/", dataset_split[2], "_no_error", 
            "/simulated_trajectory_", obs_index, ".Rdata"))
# contains Observations, Process_states

dim_proc <- dim(Process_states)
trajectory <- matrix(rep(NA, dim_proc[1] * (dim_proc[2] + 1)), nrow = dim_proc[1])
trajectory[, 1:2] <- Process_states[, 1:2]
trajectory[, 3] <- Observations[, 2]

load(paste0("data/", dataset_split[1], "/", dataset_split[2], "_with_error", 
            "/simulated_trajectory_", obs_index, ".Rdata"))
# contains Observations, pars

trajectory[, 4] <- Observations[, 2]

colnames(trajectory) <- c("time", "mRNA", "FI_no_error", "FI_with_error")


if(save_plots){
  pdf("figures_and_tables/Fig_indiv_simulated_trajectory.pdf", height = 3.2, width = 8)
}
layout(matrix(c(1,2), ncol = 2), widths = c(1,1))
par(mar=c(3.5, 4.5, 1, 1))

plot_color <-  c('#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837', '#8CB423')

plot(trajectory[,"time"], trajectory[,"mRNA"], type = 'p', pch = 20, cex = .3, #type ='l', #
     xlab = "", ylab="", las = 1, col = Helmholtz_blue)
title(ylab = bquote(X[1] ~ " (mRNA)"), cex.lab = 1.2, line = 2.8)
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)

plot(trajectory[,"time"], trajectory[,"FI_no_error"], 
     type = 'p', pch = 20, cex = .35,
     xlab = "", ylab="", las = 1, col = Helmholtz_blue)
lines(trajectory[,"time"], trajectory[,"FI_with_error"], 
      type = 'p', pch = 20, cex = .25,
      col = Helmholtz_green, lty = 1)
axis(2, labels = FALSE)
title(ylab = bquote("FI = scale * " ~ X[2] ~ " + offset"), cex.lab = 1.2, line = 2.8)
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)

legend("bottomright", pch = c(20, 20), col = c(Helmholtz_blue, Helmholtz_green),
       legend = c("data without error", "data with error"), bty ="n")


if(save_plots)dev.off()