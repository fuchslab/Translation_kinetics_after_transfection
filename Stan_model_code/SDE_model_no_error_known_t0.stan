//
// This Stan program defines a basic model of the translation kinetics
// after mRNA transfection based on an SDE model
// with known initial time point of mRNA release and
// assuming that there is no measurement error.
//

functions {
  vector drift(vector theta, vector X, real scale, real offset) {
    vector[2] drift;
    drift = [- theta[1] * X[1], 
            scale * theta[2] * X[1] - theta[3] * (X[2] - offset)]';
    return drift;
  }
  matrix diff_mat(vector theta, vector X, real scale, real offset) {
    matrix[2,2] diffusion;
    diffusion = [ [theta[1] * X[1], 0],
                        [0        , scale * (scale * theta[2] * X[1] + 
                                              theta[3] * (X[2] - offset))]];
    return diffusion;
  }
}

data {
  int<lower=0> M; // number of points observed
  real<lower=0> y_obs[M]; // observed (noisy) fluorescence intensity
  real time_points[M]; // time points of the observations
  int<lower=1, upper=M> index_before_t0; // index of the latest observed time 
  // point before the initial time point t0 when the mRNA is released
  real<lower=0> t0; // make sure that t0 is not equal to an observed time point
}

transformed data{
  // in the absence of measurement error, the offset can be directly determined 
  real<lower=0> offset; // background fluorescence 
  offset = y_obs[1];
}

parameters {
  vector<lower=0>[3] theta; // kinetic parameters (degradation rate constant of
  // mRNA, translation rate constant, degradation rate constant of GFP)
  real<lower=0> m0; // initial amount of mRNA
  real<lower=0, upper=30> scale; // scaling factor for the fluorescence signal
  real<lower=0> x_lat[M - index_before_t0]; // latent states of the diffusion process
}

transformed parameters{
  vector<lower=0>[2] x[M]; // 1st component: mRNA, 2nd component: GFP
  vector<lower=0>[2] x_t0; // initial state of the diffusion process (at t0)
  for (t in 1:(index_before_t0)) {
    x[t, 1] = 0; // before t0 there is no mRNA
  }
  x_t0 = [m0, offset]'; // at time t0 the initial amount m0 of mRNA is released
  // after t0 both components evolve as a diffusion process:
  x[(index_before_t0 + 1) : M, 1] = x_lat; 
  x[, 2] = y_obs; 
}

model {
  // factors for likelihood based on the transition probabilities of the SDE model;
  // approx. transition probability between initial time point t0 and the first 
  // observed time point after t0
  x[index_before_t0 + 1] ~
    multi_normal(x_t0 + drift(theta, x_t0, scale, offset) * (time_points[index_before_t0 + 1] - t0),
                diff_mat(theta, x_t0, scale, offset) * (time_points[index_before_t0 + 1] - t0));
  // approx. transition probabilites between observed time points
  for (t in (index_before_t0 + 2):M) {
    x[t] ~ multi_normal(x[t-1] + drift(theta, x[t-1], scale, offset) * (time_points[t] - time_points[t-1]),
                        diff_mat(theta, x[t-1], scale, offset) * (time_points[t] - time_points[t-1]));
  }
  // priors
  theta ~ normal(0, 5);
  m0 ~ normal(300,300);
}

generated quantities {
  real<lower=0> prod_theta2_m0;
  real<lower=0> prod_theta2_scale;
  real<lower=0> prod_m0_scale;
  real<lower=0> prod_theta2_m0_scale;
  prod_theta2_m0 = theta[2] * m0;
  prod_theta2_scale = theta[2] * scale;
  prod_m0_scale = m0 * scale;
  prod_theta2_m0_scale = theta[2] * m0 * scale;
}
