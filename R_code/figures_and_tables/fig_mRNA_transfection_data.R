library(readxl)
eGFP2 <- as.matrix(read_excel("data/experimental_data/eGFP/20160427_mean_eGFP.xlsx", col_names = FALSE))
d2eGFP2 <- as.matrix(read_excel("data/experimental_data/d2eGFP/20160427_mean_d2eGFP.xlsx", col_names = FALSE))


save_plots <- TRUE # TRUE FALSE

if(save_plots){
  pdf("figures_and_tables/Fig_mRNA_transfection_data.pdf", height = 3.2, width = 8)
}
layout(matrix(c(1,2), ncol = 2))
par(mar=c(3.2, 4.5, 1.5, 1))

# max of plot range for y axis
ymax <- 7500

# ---- eGFP --------------------------------------------------------------------------
# number and indices of trajectories to be ploted (>=2)
N <- 7
n <- 4
ind_eGFP <- n + (1:N)

#plot_color <- rainbow(N)
#plot_color <- c('#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')
plot_color <- c('#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')

plot(eGFP2[,1]/3600, eGFP2[,ind_eGFP[1]],  ylim = c(0, ymax), las = 1,
     type = 'p', pch = 20, cex = .3, col = plot_color[1],
     xlab = "", ylab = "")

title(main = "eGFP")
title(xlab = "time [h]", line = 2.1)
title(ylab = "fluorescence intensity [a.u.]", line = 3.5)

for(l in 2:N){
  lines(eGFP2[,1]/3600, eGFP2[ ,ind_eGFP[l]],
        type = 'p', pch = 20, cex = .3, col = plot_color[l %% length(plot_color) + 1])
}

# ---- d2eGFP ------------------------------------------------------------------------
# number and indices of trajectories to be ploted (>=2)
N <- 7
n <- 4
ind_d2eGFP <- n + (1:N)

plot(d2eGFP2[,1]/3600, d2eGFP2[ ,ind_d2eGFP[1]], ylim = c(0, ymax), las = 1,
     type = 'p', pch = 20, cex = .3, col = plot_color[1],
     xlab = "", ylab = "")

title("d2eGFP")
title(xlab = "time [h]", line = 2.1)
title(ylab = "fluorescence intensity [a.u.]", line = 3.5)

for(l in 2:N){
  lines(d2eGFP2[,1]/3600, d2eGFP2[ ,ind_d2eGFP[l]],
        type = 'p', pch = 20, cex = .3, col = plot_color[l %% length(plot_color) + 1])
}

if(save_plots) dev.off()



