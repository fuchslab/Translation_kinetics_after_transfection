library(rstan)
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

save_plots <- TRUE # TRUE FALSE

fig_height = 3
fig_width = 8.5

HMGU_red <- "#E50032"

#--- load Stanfit objects ------------------------------------------------------------
lazyLoad("preliminary_analysis/assess_need_for_data_augmentation/mRNA1_no_error_no_scale_with_imputation_latent_var_different_m_cache/html/stanfit1i_e6ce07fb503b98f3a3688dedd155ed71")
lazyLoad("preliminary_analysis/assess_need_for_data_augmentation/mRNA1_no_error_no_scale_with_imputation_latent_var_different_m_cache/html/stanfit2i_a5958ec413dfa072a09c17ceb2f6aa9d")
lazyLoad("preliminary_analysis/assess_need_for_data_augmentation/mRNA1_no_error_no_scale_with_imputation_latent_var_different_m_cache/html/stanfit5i_49ca3afcc41667b4df9ef4b12d3d6587")
lazyLoad("preliminary_analysis/assess_need_for_data_augmentation/mRNA1_no_error_no_scale_with_imputation_latent_var_different_m_cache/html/stanfit10i_4e0b40a8ba2c02da9d4ad39ab7b74088")

m_vec <- c(1,2,5,10)
calc_CI <- function(stanfit_object){
  res_mat <- matrix(NA, nrow = 3, ncol = 3)
  rownames(res_mat) <- c("theta[1]", "theta[2]", "theta[3]")
  colnames(res_mat) <- c("median", "hpdi_l", "hpdi_u")
  theta <- rstan::extract(stanfit_object, pars="theta")[[1]]
  res_mat[1,] <- quantile(theta[,1], probs = c(.5, .025, .975), type = 7)
  res_mat[2,] <- quantile(theta[,2], probs = c(.5, .025, .975), type = 7)
  res_mat[3,] <- quantile(theta[,3], probs = c(.5, .025, .975), type = 7)
  return(res_mat)
}

HMC_m1 <- calc_CI(mRNA1_no_error_no_scale_latent_var_provide_init_no_imputation)
HMC_m2 <- calc_CI(mRNA1_no_error_no_scale_with_imputation_init_val2)
HMC_m5 <- calc_CI(mRNA1_no_error_no_scale_with_imputation_init_val5)
HMC_m10 <- calc_CI(mRNA1_no_error_no_scale_with_imputation_init_val10)

parameter <- rep(c("theta[1]", "theta[2]", "theta[3]"), times = length(m_vec))
m_value <- rep(as.factor(m_vec), each = 3)

res <- rbind(HMC_m1, HMC_m2, HMC_m5, HMC_m10)
df <- data.frame(m_value, parameter, res)

#-- genereate the plot ---------------------------------------------------------------
if(save_plots){
  pdf("figures_and_tables/Fig_assess_need_for_data_augm_with_Stan.pdf", width = fig_width, height = fig_height)
}

res_theta1 <- df %>% filter(parameter == "theta[1]")
plot_theta1 <- ggplot(data=res_theta1, aes(x = m_value, y = median)) +
  geom_errorbar(aes(ymax = hpdi_l, ymin = hpdi_u), width=.5, color = HMGU_red) +
  geom_point(aes(y=median), shape=20, size=3, color = HMGU_red) +
  geom_abline(intercept=.11, slope=0, color="black", lwd=.4) +
  ylab("Point estimate and 95% CI") +
  ggtitle("theta[1]") +
  xlab("")+ theme(legend.position="none",
                  plot.title = element_text(hjust = .5))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=13),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=12))

res_theta2 <- df %>% filter(parameter == "theta[2]")
plot_theta2 <- ggplot(data=res_theta2, aes(x = m_value, y = median)) +
  geom_errorbar(aes(ymax = hpdi_l, ymin = hpdi_u), width=.5, color = HMGU_red) +
  geom_point(aes(y=median), shape=20, size=3, color = HMGU_red) +
  geom_abline(intercept=.3, slope=0, color="black", lwd=.4) +
  ylab("") +
  ggtitle("theta[2]") +
  xlab("# of inter-obs. intervals") +
  theme(legend.position="none", plot.title = element_text(hjust = .5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=13),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=12),
        panel.grid.major.x = element_blank())

res_theta3 <- df %>% filter(parameter == "theta[3]")
plot_theta3 <- ggplot(data=res_theta3, aes(x = m_value, y = median)) +
  geom_errorbar(aes(ymax = hpdi_l, ymin = hpdi_u), width=.5, color = HMGU_red) +
  geom_point(aes(y=median), shape=20, size=3, color = HMGU_red) +
  geom_abline(intercept=.09, slope=0, color="black", lwd=.4) +
  ggtitle("theta[3]") +
  ylab("") +
  xlab("")+ theme(legend.position="none",
                  plot.title = element_text(hjust = .5))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=13),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=12))
gridExtra::grid.arrange(plot_theta1, plot_theta2, plot_theta3, ncol = 3,
                        widths = c( unit(5, c("cm")), unit(5, c("cm")), unit(5, c("cm"))))

if(save_plots) dev.off()

# extract the computational times ----------------------------------------------------
# time_m1 <- get_elapsed_time(mRNA1_no_error_no_scale_latent_var_provide_init_no_imputation)
# time_m2 <- get_elapsed_time(mRNA1_no_error_no_scale_with_imputation_init_val2)
# time_m5 <- get_elapsed_time(mRNA1_no_error_no_scale_with_imputation_init_val5)
# time_m10 <-get_elapsed_time(mRNA1_no_error_no_scale_with_imputation_init_val10)
#
# comp_times_stan <- matrix(NA, ncol = 2, nrow = length(m_vec))
# rownames(comp_times_stan) <- paste("m = ", m_vec)
# colnames(comp_times_stan) <- c("mean", "sd")
# comp_times_stan[1,] <- c(mean(rowSums(time_m1)), sd(rowSums(time_m1)))
# comp_times_stan[2,] <- c(mean(rowSums(time_m2)), sd(rowSums(time_m2)))
# comp_times_stan[3,] <- c(mean(rowSums(time_m5)), sd(rowSums(time_m5)))
# comp_times_stan[4,] <- c(mean(rowSums(time_m10)), sd(rowSums(time_m10)))
# comp_times_stan

