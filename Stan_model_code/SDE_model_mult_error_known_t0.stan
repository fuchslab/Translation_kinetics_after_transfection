//
// This Stan program defines a basic model of the translation kinetics
// after mRNA transfection based on an SDE model
// with known initial time point of mRNA release and
// assuming multiplicative normal measurement errors.
//

functions {
  vector drift(vector theta, vector X) {
    vector[2] drift;
    drift = [- theta[1] * X[1], theta[2] * X[1] - theta[3] * X[2] ]';
    return drift;
  }
  matrix diff_mat(vector theta, vector X) {
    matrix[2,2] diffusion;
    diffusion = [ [theta[1] * X[1], 0],
                        [0        , theta[2] * X[1] + theta[3] * X[2]]];
    return diffusion;
  }
}

data {
  int<lower=0> M; // number of points observed
  vector<lower=0>[M] y_obs; // observed (noisy) fluorescence intensity
  real time_points[M]; // time points of the observations
  int<lower=1, upper=M> index_before_t0; // index of the latest observed time 
  // point before the initial time point t0 when the mRNA is released
  real<lower=0> t0; // make sure that t0 is not equal to an observed time point
}

transformed data{
  // log-transform the observed data, as we assume multiplicative measurement errors
  vector[M] log_y_obs;
  log_y_obs = log(y_obs);
}

parameters {
  vector<lower=0>[3] theta; // kinetic parameters (degradation rate constant of
  // mRNA, translation rate constant, degradation rate constant of GFP)
  real<lower=0> m0; // initial amount of mRNA
  real<lower=1e-3, upper=10> sigma; // standard deviation of measurement error
  real<lower=0, upper=30> scale; // scaling factor for the fluorescence signal
  real<lower=0, upper=30> offset; // background fluorescence
  vector<lower=0>[2] x_k[M - index_before_t0]; // latent states of the diffusion process
}

transformed parameters{
  vector<lower=0>[2] x[M]; // 1st component: mRNA, 2nd component: GFP
  vector<lower=0>[2] x_t0; // initial state of the diffusion process (at t0)
  for (t in 1:(index_before_t0)) {
    x[t] = [0,0]'; // before t0 both components are zero
  }
  x_t0 = [m0,0]'; // at time t0 the initial amount m0 of mRNA is released
  x[(index_before_t0 + 1) : M] = x_k; // after t0 both components evolve as 
  // a diffusion process
}

model {
  // factors for likelihood based on the transition probabilities of the SDE model;
  //approx. transition probability between initial time point t0 and the first 
  // observed time point after t0
  x[index_before_t0 + 1] ~
    multi_normal(x_t0 + drift(theta, x_t0) * (time_points[index_before_t0 + 1] - t0),
                diff_mat(theta, x_t0) * (time_points[index_before_t0 + 1] - t0));
  // approx. transition probabilites between observed time points
  for (t in (index_before_t0 + 2):M) {
    x[t] ~ multi_normal(x[t-1] + drift(theta, x[t-1]) * (time_points[t] - time_points[t-1]),
                        diff_mat(theta, x[t-1]) * (time_points[t] - time_points[t-1]));
  }
  // factors for likelihood assuming a linear transformation of the number of protein
  // molecules into the fluorescence signal and
  // multiplicative normal measurement errors
  for (i in 1:M){
    log_y_obs[i] ~ normal(log(scale * x[i,2] + offset), sigma);
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
