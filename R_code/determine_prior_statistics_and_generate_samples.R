# This R script is used to calculate the 95% interquantile range of the prior 
# distributions used for each parameter

param <- c("theta[1]", "theta[2]", "theta[3]", "m0", "scale", "offset", "t0",
           "sigma", "prod_theta2_m0", 
           "prod_theta2_scale", "prod_m0_scale", "prod_theta2_m0_scale")
col_names <- c("mean", "2.5%", "50%", "97.5%", "length of 95% HPDI")

n_col <- length(col_names)
n_row <- length(param)

res_mat <- matrix(rep(NA, n_col * n_row), ncol = n_col)
colnames(res_mat) <- col_names
rownames(res_mat) <- param

probs <- c(0.025, 0.5, 0.975)

# theta ~ truncated_Normal>=0 (mean = 0, sd = 5)--------------------------------
a <- 0
mu <- 0
sig <- 5
mean_theta <- mu + dnorm(a, mean = mu, sd = sig) / 
  (1 - pnorm(a, mean = mu, sd = sig)) * sig
quantiles <- msm::qtnorm(p = probs, mean = mu, sd = sig, lower = a)
len_95_hpdi <- quantiles[3] - quantiles[1]
 
res_mat["theta[1]", ] <-
  res_mat["theta[2]", ] <-
  res_mat["theta[3]", ] <- c(mean_theta, quantiles, len_95_hpdi)

# m0 ~ truncated_Normal>=0 (mean = 300, sd = 300)--------------------------------
a <- 0
mu <- 300
sig <- 300
mean <- mu + dnorm(a, mean = mu, sd = sig) / 
  (1 - pnorm(a, mean = mu, sd = sig)) * sig
quantiles <- msm::qtnorm(p = probs, mean = mu, sd = sig, lower = a)
len_95_hpdi <- quantiles[3] - quantiles[1]

res_mat["m0", ] <- c(mean, quantiles, len_95_hpdi)

# scale, offset, t0 ~ Unif(0, 30)-----------------------------------------------
a <- 0
b <- 30
mean <- (a + b) / 2
quantiles <- qunif(p = probs, min = a, max = b)
len_95_hpdi <- quantiles[3] - quantiles[1]

res_mat["scale", ] <- 
  res_mat["offset", ] <- 
  res_mat["t0", ] <- c(mean, quantiles, len_95_hpdi)

# sigma ~ Unif(1e-3, 10)--------------------------------------------------------
a <- 1e-3
b <- 10
mean <- (a + b) / 2
quantiles <- qunif(p = probs, min = a, max = b)
len_95_hpdi <- quantiles[3] - quantiles[1]

res_mat["sigma", ] <- c(mean, quantiles, len_95_hpdi)

## products --------------------------------------------------------------------
set.seed(1569791933)
N <- 1e8

sample_theta_2 <- msm::rtnorm(N, mean = 0, sd = 5, lower = 0)
sample_m0 <- msm::rtnorm(N, mean = 300, sd = 300, lower = 0)
sample_scale <- runif(N, min = 0, max = 30)

sample_prod_theta2_m0 <- sample_theta_2 * sample_m0
mean <- mean(sample_prod_theta2_m0)
quantiles <- quantile(sample_prod_theta2_m0, p = probs)
len_95_hpdi <- quantiles[3] - quantiles[1]
res_mat["prod_theta2_m0", ] <- c(mean, quantiles, len_95_hpdi)

sample_prod_theta2_scale <- sample_theta_2 * sample_scale
mean <- mean(sample_prod_theta2_scale)
quantiles <- quantile(sample_prod_theta2_scale, p = probs)
len_95_hpdi <- quantiles[3] - quantiles[1]
res_mat["prod_theta2_scale", ] <- c(mean, quantiles, len_95_hpdi)

sample_prod_m0_scale <- sample_m0 * sample_scale
mean <- mean(sample_prod_m0_scale)
quantiles <- quantile(sample_prod_m0_scale, p = probs)
len_95_hpdi <- quantiles[3] - quantiles[1]
res_mat["prod_m0_scale", ] <- c(mean, quantiles, len_95_hpdi)

sample_prod_theta2_m0_scale <- sample_theta_2 * sample_m0 * sample_scale
mean <- mean(sample_prod_theta2_m0_scale)
quantiles <- quantile(sample_prod_theta2_m0_scale, p = probs)
len_95_hpdi <- quantiles[3] - quantiles[1]
res_mat["prod_theta2_m0_scale", ] <- c(mean, quantiles, len_95_hpdi)

saveRDS(res_mat, 
        file = "intermediate_output_files/prior_statistics.rds")

# construct list of prior density functions and samples ------------------------
# for parameters where density function of the prior distribution is unknown
# save (part of) the sample generated before
N_2 <- 1e7
prior_dens_or_sample <-
  list("theta[1]" = list(dens_fct_exists = TRUE,
                         densfct = function(x){
                           msm::dtnorm(x, mean = 0, sd = 5, lower = 0)}),
       "theta[2]" = list(dens_fct_exists = TRUE,
                         densfct = function(x){
                           msm::dtnorm(x, mean = 0, sd = 5, lower = 0)}),
       "theta[3]" = list(dens_fct_exists = TRUE,
                         densfct = function(x){
                           msm::dtnorm(x, mean = 0, sd = 5, lower = 0)}),
       "m0" = list(dens_fct_exists = TRUE,
                         densfct = function(x){
                           msm::dtnorm(x, mean = 300, sd = 300, lower = 0)}),
       "scale" = list(dens_fct_exists = TRUE,
                   densfct = function(x){ dunif(x, min = 0, max = 30)}),
       "offset" = list(dens_fct_exists = TRUE,
                      densfct = function(x){ dunif(x, min = 0, max = 30)}),
       "t0" = list(dens_fct_exists = TRUE,
                      densfct = function(x){ dunif(x, min = 0, max = 30)}),
       "sigma" = list(dens_fct_exists = TRUE,
                   densfct = function(x){ dunif(x, min = 1e-3, max = 10)}),
       "prod_theta2_m0" = list(dens_fct_exists = FALSE,
                               sample = sample_prod_theta2_m0[1:N_2]),
       "prod_theta2_scale" = list(dens_fct_exists = FALSE,
                                  sample = sample_prod_theta2_scale[1:N_2]),
       "prod_m0_scale" = list(dens_fct_exists = FALSE,
                              sample = sample_prod_m0_scale[1:N_2]),
       "prod_theta2_m0_scale" = list(dens_fct_exists = FALSE,
                                     sample = sample_prod_theta2_m0_scale[1:N_2]))

saveRDS(prior_dens_or_sample, 
        file = "intermediate_output_files/prior_dens_fct_or_sample.rds")