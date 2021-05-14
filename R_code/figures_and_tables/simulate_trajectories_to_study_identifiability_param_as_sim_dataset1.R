# This R script is used to simulate trajectories of the SDE and the ODE model
# in order to study parameter identifiability.

source("R_code/functions_to_simulate_dynamical_systems/create_model.R")
source("R_code/functions_to_simulate_dynamical_systems/diffusion.R")

save_plots <- TRUE #TRUE #FALSE

# ODE solution -----------------------------------------------------------------
ODE_solution <- function(time, initial, theta, t0){
  ODE_sol <- time
  ODE_sol[time < t0] <- 0
  indizes_after_t0 <- time > t0
  ODE_sol[indizes_after_t0] <- theta[2] * initial[1] / (theta[3] - theta[1]) *
    (exp(-theta[1] * (time[indizes_after_t0] - t0)) - exp(-theta[3] * (time[indizes_after_t0] - t0)))
  ODE_sol
}


# mRNA transfection - model 1 --------------------------------------------------

num_paths <- 30 # number of trajectories to generate

delta_t <- .01
Time <- 30
time <- seq(from = 0, to = Time, by = delta_t)

t0 <- 0 #0.96
offset <- 6.5


# prod_theta2_m0_scale: only the product of these three parameters is
# identifiable in the ODE model and therefore is kept constant here
prod_theta2_m0_scale <- 241.92


##----- parameter setting 1 -----------------------------------------------------------
theta2_1 <- 0.32 # translation rate constant 
m0_1 <- 420 # initial amount of mRNA
scale_1 <- prod_theta2_m0_scale / (theta2_1 * m0_1) # scaling factor for the fluorescence

# initial conditions for mRNA, GFP
initial <- c(m0_1, 0)
# parameters:
pars1 <- list(theta = c(0.2, theta2_1, 0.01), m0 = m0_1, t0 = t0, 
                      scale = scale_1, offset = offset)

# mRNA transfection - model 1
# initial conditions for mRNA, GFP
initial <- c(pars1$m0, 0)
# parameters:
kinetic_pars <- pars1$theta
# now create a model object:
model1 <- create_model_mRNA1(initial = initial, theta = kinetic_pars)

ODE1 <- ODE_solution(time, initial, kinetic_pars, pars1$t0)

##----- parameter setting 2 -----------------------------------------------------------
theta2_2 <- theta2_1 / 10 # translation rate costant
m0_2 <- m0_1 # initial amount of mRNA
scale_2 <- prod_theta2_m0_scale / (theta2_2 * m0_2) # scaling factor for the fluorescence

# parameters:
pars2 <- list(theta = c(0.2, theta2_2, 0.01), m0 = m0_2, t0 = t0, 
              scale = scale_2, offset = offset)

# initial conditions for mRNA, GFP
initial <- c(pars2$m0, 0)
# parameters:
kinetic_pars <- pars2$theta
# now create a model object:
model2 <- create_model_mRNA1(initial = initial, theta = kinetic_pars)

ODE2 <- ODE_solution(time, initial, kinetic_pars, pars2$t0)

##----- parameter setting 3 -----------------------------------------------------------
theta2_3 <- theta2_1 # translation rate constant
m0_3 <- m0_1 * 10 # initial amount of mRNA
scale_3 <- prod_theta2_m0_scale / (theta2_3 * m0_3) # scaling factor for the fluorescence

# parameters:
pars3 <- list(theta = c(0.2, theta2_3, 0.01), m0 = m0_3, t0 = t0, 
              scale = scale_3, offset = offset)

# initial conditions for mRNA, GFP
initial <- c(pars3$m0, 0)
# parameters:
kinetic_pars <- pars3$theta
# now create a model object:
model3 <- create_model_mRNA1(initial = initial, theta = kinetic_pars)

ODE3 <- ODE_solution(time, initial, kinetic_pars, pars3$t0)

##----- parameter setting 4 -----------------------------------------------------------
theta2_4 <- theta2_1 / 10 # translation rate constant
m0_4 <-  m0_1 * 10 # initial amount of mRNA
scale_4 <- prod_theta2_m0_scale / (theta2_4 * m0_4) # scaling factor for the fluorescence

# parameters:
pars4 <- list(theta = c(0.2, theta2_4, 0.01), m0 = m0_4, t0 = t0, 
              scale = scale_4, offset = offset)

# initial conditions for mRNA, GFP
initial <- c(pars4$m0, 0)
# parameters:
kinetic_pars <- pars4$theta
# now create a model object:
model4 <- create_model_mRNA1(initial = initial, theta = kinetic_pars)

ODE4 <- ODE_solution(time, initial, kinetic_pars, pars4$t0)

##----- parameter setting 5 -----------------------------------------------------------
theta2_5 <- theta2_1 * 100 # translation rate constant
m0_5 <- m0_1  # initial amount of mRNA
scale_5 <- prod_theta2_m0_scale / (theta2_5 * m0_5) # scaling factor for the fluorescence

# parameters:
pars5 <- list(theta = c(0.2, theta2_5, 0.01), m0 = m0_5, t0 = t0, 
              scale = scale_5, offset = offset)

# initial conditions for mRNA, GFP
initial <- c(pars5$m0, 0)
# parameters:
kinetic_pars <- pars5$theta
# now create a model object:
model5 <- create_model_mRNA1(initial = initial, theta = kinetic_pars)

ODE5 <- ODE_solution(time, initial, kinetic_pars, pars5$t0)

##----- parameter setting 6 -----------------------------------------------------------
theta2_6 <- theta2_1 * 10 # translation rate constant
m0_6 <- m0_1 / 10 # initial amount of mRNA
scale_6 <- prod_theta2_m0_scale / (theta2_6 * m0_6) # scaling factor for the fluorescence

# parameters:
pars6 <- list(theta = c(0.2, theta2_6, 0.01), m0 = m0_6, t0 = t0, 
              scale = scale_6, offset = offset)

# initial conditions for mRNA, GFP
initial <- c(pars6$m0, 0)
# parameters:
kinetic_pars <- pars6$theta
# now create a model object:
model6 <- create_model_mRNA1(initial = initial, theta = kinetic_pars)

ODE6 <- ODE_solution(time, initial, kinetic_pars, pars6$t0)

##-- simulation keeping the product theta2*m0*scale constant -------------------
mRNA1 <- matrix(NA, ncol = length(time), nrow = num_paths)
GFP1 <- matrix(NA, ncol = length(time), nrow = num_paths)
mRNA2 <- matrix(NA, ncol = length(time), nrow = num_paths)
GFP2 <- matrix(NA, ncol = length(time), nrow = num_paths)
mRNA3 <- matrix(NA, ncol = length(time), nrow = num_paths)
GFP3 <- matrix(NA, ncol = length(time), nrow = num_paths)
mRNA4 <- matrix(NA, ncol = length(time), nrow = num_paths)
GFP4 <- matrix(NA, ncol = length(time), nrow = num_paths)
mRNA5 <- matrix(NA, ncol = length(time), nrow = num_paths)
GFP5 <- matrix(NA, ncol = length(time), nrow = num_paths)
mRNA6 <- matrix(NA, ncol = length(time), nrow = num_paths)
GFP6 <- matrix(NA, ncol = length(time), nrow = num_paths)


seed <- 1942308132
set.seed(seed)
for(i in 1:num_paths){
  process <- diffusion(model1, maxtime=30, tstep = delta_t)
  mRNA1[i,] <- process[,2]
  GFP1[i,] <- process[,3]
}

set.seed(seed)
for(i in 1:num_paths){
  process <- diffusion(model2, maxtime=30, tstep = delta_t)
  mRNA2[i,] <- process[,2]
  GFP2[i,] <- process[,3]
}

set.seed(seed)
for(i in 1:num_paths){
  process <- diffusion(model3, maxtime=30, tstep = delta_t)
  mRNA3[i,] <- process[,2]
  GFP3[i,] <- process[,3]
}

set.seed(seed)
for(i in 1:num_paths){
  process <- diffusion(model4, maxtime=30, tstep = delta_t)
  mRNA4[i,] <- process[,2]
  GFP4[i,] <- process[,3]
}

set.seed(seed)
for(i in 1:num_paths){
  process <- diffusion(model5, maxtime=30, tstep = delta_t)
  mRNA5[i,] <- process[,2]
  GFP5[i,] <- process[,3]
}

set.seed(seed)
for(i in 1:num_paths){
  process <- diffusion(model6, maxtime=30, tstep = delta_t)
  mRNA6[i,] <- process[,2]
  GFP6[i,] <- process[,3]
}

FI1 <- GFP1*scale_1 + offset
FI2 <- GFP2*scale_2 + offset
FI3 <- GFP3*scale_3 + offset
FI4 <- GFP4*scale_4 + offset
FI5 <- GFP5*scale_5 + offset
FI6 <- GFP6*scale_6 + offset

max_FI <- max(FI1, FI2, FI3, FI4, FI5, FI6)

FI_ODE1 <- ODE1 * scale_1 + offset
FI_ODE2 <- ODE2 * scale_2 + offset
FI_ODE3 <- ODE3 * scale_3 + offset
FI_ODE4 <- ODE4 * scale_4 + offset
FI_ODE5 <- ODE5 * scale_5 + offset
FI_ODE6 <- ODE6 * scale_6 + offset

##-------- plot all 6 parameter settings --------------------------------------------
if(save_plots){
  pdf("figures_and_tables/Fig_study_identifiability_keep_prod_theta2_m0_scale_constant.pdf", width = 7.5, height = 4.5)
}

left_margin_no_labels <- 1.5

layout(matrix(c(1,2,3,4,5,6), ncol = 3), widths = c(1.15,1,1))
par(mar=c(3.5, 4.3, 2, 0.5))

plot_color <-  rep(c('#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837', '#8CB423'), times = ceiling(num_paths/8))
#excluded because too dark'#004529'

plot(process[,1], FI1[1,], type = 'l', xlab = "", ylab="",
     las = 1, col = plot_color[1], ylim = c(1, max_FI))
# title(ylab = "fluor. intensity = scale * GFP + offset", cex.lab = 1.2, line = 3.3)
title(ylab = bquote("scale *" ~ X[2] ~ " + offset"), cex.lab = 1.2, line = 3)
#title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)
title(main = bquote(~m[0]==.(m0_1)~", " ~ theta[2] == .(theta2_1) ~ ", scale = " ~ .(scale_1)))
for(l in 2:num_paths){
  lines(time, FI1[l,], col = plot_color[l])
}

lines(time, FI_ODE1, col = 1, lwd = 2, lty = 1)

plot(process[,1], FI2[1,], type = 'l', xlab = "", ylab="",
     las = 1, col = plot_color[1], ylim = c(1, max_FI))
# title(ylab = "fluor. intensity = scale * GFP + offset", cex.lab = 1.2, line = 3.3)
title(ylab = bquote("scale *" ~ X[2] ~ " + offset"), cex.lab = 1.2, line = 3)
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)
title(main = bquote(~ m[0] == .(m0_2) ~ ", " ~ theta[2] == .(theta2_2) ~ ", scale = " ~ .(scale_2)))
for(l in 2:num_paths){
  lines(time, FI2[l,], col = plot_color[l])
}

lines(time, FI_ODE2, col = 1, lwd = 2, lty = 1)


par(mar=c(3.5, left_margin_no_labels, 2, 0.5))
plot(process[,1], FI3[1,], type = 'l', xlab = "", ylab="",  yaxt='n',
     las = 1, col = plot_color[1], ylim = c(1, max_FI))
axis(2, labels = FALSE)
#title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)
title(main = bquote(~ m[0] == .(m0_3) ~ ", " ~ theta[2] == .(theta2_3) ~ ", scale = " ~ .(scale_3)))
for(l in 2:num_paths){
  lines(time, FI3[l,], col = plot_color[l])
}

lines(time, FI_ODE3, col = 1, lwd = 2, lty = 1)

plot(process[,1], FI4[1,], type = 'l', xlab = "", ylab="",  yaxt='n',
     las = 1, col = plot_color[1], ylim = c(1, max_FI))
axis(2, labels = FALSE)
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)
title(main = bquote(~ m[0] == .(m0_4) ~ ", " ~ theta[2] == .(theta2_4) ~ ", scale = " ~ .(scale_4)))
for(l in 2:num_paths){
  lines(time, FI4[l,], col = plot_color[l])
}

lines(time, FI_ODE4, col = 1, lwd = 2, lty = 1)


plot(process[,1], FI5[1,], type = 'l', xlab = "", ylab="",  yaxt='n',
     las = 1, col = plot_color[1], ylim = c(1, max_FI))
axis(2, labels = FALSE)
#title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)
title(main = bquote(~ m[0] == .(m0_5) ~ ", " ~ theta[2] == .(theta2_5) ~ ", scale = " ~ .(scale_5)))
for(l in 2:num_paths){
  lines(time, FI5[l,], col = plot_color[l])
}

lines(time, FI_ODE5, col = 1, lwd = 2, lty = 1)

h <- 2
v <- 320
legend(x = 10+h, y = 35+v, bty = "n", legend = c("ODE solution", "SDE trajectories"),
       col = c(1, plot_color[3]), lty = c(1,1), lwd = c(2,1))
legend(x = 10+h, y = 65+v, bty = "n", legend = c("", ""), col = c(0, plot_color[2]), lty = c(1,1))
legend(x = 10+h, y = 5+v, bty = "n", legend = c("", ""), col = c(0, plot_color[7]), lty = c(1,1))
legend(x = 10+h, y = -25+v, bty = "n", legend = c("", ""), col = c(0, plot_color[6]), lty = c(1,1))


plot(process[,1], FI6[1,], type = 'l', xlab = "", ylab="",  yaxt='n',
     las = 1, col = plot_color[1], ylim = c(1, max_FI))
axis(2, labels = FALSE)
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)
title(main = bquote(~ m[0] == .(m0_6) ~ ", " ~ theta[2] == .(theta2_6) ~ ", scale = " ~ .(scale_6)))
for(l in 2:num_paths){
  lines(time, FI6[l,], col = plot_color[l])
}

lines(time, FI_ODE6, col = 1, lwd = 2, lty = 1)

if(save_plots) dev.off()


##----- parameter setting 7 -----------------------------------------------------------
# parameters:
pars7 <- list(theta = c(0.01, theta2_1, 0.2), m0 = m0_1, t0 = t0, 
              scale = scale_1, offset = offset)

# initial conditions for mRNA, GFP
initial <- c(pars7$m0, 0)
# parameters:
kinetic_pars <- pars7$theta
# now create a model object:
model7 <- create_model_mRNA1(initial = initial, theta = kinetic_pars)

ODE7 <- ODE_solution(time, initial, kinetic_pars, pars7$t0)
##-------- simulation swapping theta1 and theta3-------------------------------
mRNA7 <- matrix(NA, ncol = length(time), nrow = num_paths)
GFP7 <- matrix(NA, ncol = length(time), nrow = num_paths)

set.seed(seed)
for(i in 1:num_paths){
  process <- diffusion(model7, maxtime=30, tstep = delta_t)
  mRNA7[i,] <- process[,2]
  GFP7[i,] <- process[,3]
}

FI7 <- GFP7*scale_1 + offset

max_FI2 <- max(FI1, FI7)


FI_ODE7 <- ODE7 * scale_1 + offset

##-------- plot 3 -- switch theta_1 and theta_3 --------------------------------------------------------------------
if(save_plots){
  pdf("figures_and_tables/Fig_study_identifiability_swap_theta1_theta3_several_trajectories.pdf", width = 7, height = 3)
}

layout(matrix(c(1,2), ncol = 2), widths = c(1.15,1,1))
par(mar=c(3.5, 4, 2, 1))

plot_color <-  rep(c('#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837', '#8CB423'), times = ceiling(num_paths/8))

plot(process[,1], FI1[1,], type = 'l', xlab = "", ylab="",
     las = 1, col = plot_color[1], ylim = c(1, max_FI2))
# title(ylab = "fluor. intensity = scale * GFP + offset", cex.lab = 1.2, line = 3)
title(ylab = bquote("scale * " ~ X[2] ~ " + offset"), cex.lab = 1.2, line = 2.8)
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)
title(main = bquote(~ theta[1] == .(pars1$theta[1]) ~ ", "  ~ theta[3] == .(pars1$theta[3])))
for(l in 2:num_paths){
  lines(time, FI1[l,], col = plot_color[l])
}

lines(time, FI_ODE1, col = 1, lwd = 2, lty = 1)

par(mar=c(3.5, 1.5, 2, 1))
plot(process[,1], FI7[1,], type = 'l', xlab = "", ylab="",  yaxt='n',
     las = 1, col = plot_color[1], ylim = c(1, max_FI2))
axis(2, labels = FALSE)
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)
title(main = bquote(~ theta[1] == .(pars7$theta[1]) ~ ", "  ~ theta[3] == .(pars7$theta[3])))
for(l in 2:num_paths){
  lines(time, FI7[l,], col = plot_color[l])
}

lines(time, FI_ODE7, col = 1, lwd = 2, lty = 1)

h <- -1
v <- 320
legend(x = 10+h, y = 40+v, bty = "n", legend = c("ODE solution", "SDE trajectories"),
       col = c(1, plot_color[3]), lty = c(1,1), lwd = c(2,1))
legend(x = 10+h, y = 60+v, bty = "n", legend = c("", ""), col = c(0, plot_color[2]), lty = c(1,1))
legend(x = 10+h, y = 20+v, bty = "n", legend = c("", ""), col = c(0, plot_color[7]), lty = c(1,1))
legend(x = 10+h, y = 0+v, bty = "n", legend = c("", ""), col = c(0, plot_color[6]), lty = c(1,1))

if(save_plots)dev.off()

##-------- plot 4 -- switch theta_1 and theta_3 --------------------------------------------------------------------
if(save_plots){
  pdf("figures_and_tables/Fig_study_identifiability_swap_theta1_theta3_individual_trajectories.pdf", width = 7.5, height = 4.25)
}

layout(matrix(c(1,2,3,4,5,6), ncol = 3), widths = c(1.15,1,1))
par(mar=c(3.5, 4.3, 1, 1))

plot_color <-  rep(c('#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529',  '#8CB423'), times = ceiling(num_paths/8))

plot(process[,1], FI1[1,], type = 'l', xlab = "", ylab="",
     las = 1, col = plot_color[2], ylim = c(1, max_FI2))
lines(process[,1], FI7[1,], type = 'l', col = plot_color[6])
# title(ylab = "fluor. intensity = scale * GFP + offset", cex.lab = 1.2, line = 3.3)
title(ylab = bquote("scale * " ~ X[2] ~ " + offset"), cex.lab = 1.2, line = 3)

plot(process[,1], FI1[2,], type = 'l', xlab = "", ylab="",
     las = 1, col = plot_color[2], ylim = c(1, max_FI2))
lines(process[,1], FI7[2,], type = 'l', col = plot_color[6])
# title(ylab = "fluor. intensity = scale * GFP + offset", cex.lab = 1.2, line = 3.3)
title(ylab = bquote("scale *" ~ X[2] ~ " + offset"), cex.lab = 1.2, line = 3)
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)


par(mar=c(3.5, 1.5, 1, 1))
plot(process[,1], FI1[3,], type = 'l', xlab = "", ylab="",  yaxt='n',
     las = 1, col = plot_color[2], ylim = c(1, max_FI2))
axis(2, labels = FALSE)
lines(process[,1], FI7[3,], type = 'l', col = plot_color[6])

plot(process[,1], FI1[4,], type = 'l', xlab = "", ylab="",  yaxt='n',
     las = 1, col = plot_color[2], ylim = c(1, max_FI2))
axis(2, labels = FALSE)
lines(process[,1], FI7[4,], type = 'l', col = plot_color[6])
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)

plot(process[,1], FI1[5,], type = 'l', xlab = "", ylab="",  yaxt='n',
     las = 1, col = plot_color[2], ylim = c(1, max_FI2))
axis(2, labels = FALSE)
lines(process[,1], FI7[5,], type = 'l', col = plot_color[6])

plot(process[,1], FI1[6,], type = 'l', xlab = "", ylab="",  yaxt='n',
     las = 1, col = plot_color[2], ylim = c(1, max_FI2))
axis(2, labels = FALSE)
lines(process[,1], FI7[6,], type = 'l', col = plot_color[6])
title(xlab = "time [h]", cex.lab = 1.2, line = 2.5)

legend("bottomright",
       legend = c(as.expression(bquote(~ theta[1] == .(pars1$theta[1]) ~ ",   "  ~ theta[3] == .(pars1$theta[3]))),
                  as.expression(bquote(~ theta[1] == .(pars7$theta[1]) ~ ", "  ~ theta[3] == .(pars7$theta[3])))),
       lty = c(1,1), lwd = c(1.5,1.5), col = c(plot_color[2], plot_color[6]),
       bty = "n", cex = 1.15)


if(save_plots) dev.off()
