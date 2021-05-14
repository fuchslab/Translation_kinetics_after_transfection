//
// This Stan program defines a basic model of the translation kinetics
// after mRNA transfection based on an ODE model with explicit solution
// with unknown initial time point of mRNA release and
// assuming multiplicative normal measurement errors.
//

functions {
  // function to obtain the amount of GFP at the given time points based on 
  // the solution of the ODE model
  vector ode_sol(real t0,
  real[] time_points,
  int M,
  real m0,
  real[] theta) {
    vector[M] x_sol; // protein (GFP)
    
    for (t in 1:M){
      if(time_points[t] < t0){
        x_sol[t] = 0;
      }else{
        if(theta[1] != theta[3]){
          x_sol[t] = theta[2] * m0 / (theta[3] - theta[1]) *
          (exp(-theta[1] * (time_points[t] - t0)) - exp(-theta[3] * (time_points[t] - t0)));
        }else{
          x_sol[t] = theta[2] * m0 *
          (time_points[t] - t0) * exp(-theta[3] * (time_points[t] - t0));
        }
      }
    }
    return x_sol;
  }
}

data {
  int<lower=0> M; // number of points observed
  vector<lower=0>[M] y_obs; // observed (noisy) fluorescence intensity
  real time_points[M]; // time points of the observations
}

transformed data{
  // log-transform the observed data, as we assume multiplicative measurement errors
  vector[M] log_y_obs;
  log_y_obs = log(y_obs);
}

parameters {
  real<lower=0> theta[3]; // kinetic parameters (degradation rate constant of 
  // mRNA, translation rate constant, degradation rate constant of GFP)
  real<lower=0> m0; // initial amount of mRNA
  real<lower=1e-3, upper=10> sigma; // standard deviation of measurement error
  real<lower=0, upper=30> scale; // scaling factor for the fluorescence signal
  real<lower=0, upper=30> offset; // background fluorescence
  real<lower=0, upper=time_points[M]> t0; // time point at which inital amount 
  // of mRNA is released
}

model {
  vector[M] log_y_sim;
  // first linear transformation of the ODE solution for the GFP to the 
  // fluorescence signal; then log transformation of the fluorescence signal
  log_y_sim = log(scale * ode_sol(t0, time_points, M, m0, theta) + 
              rep_vector(offset, M));
  // likelihood assuming multiplicative normal measurement errors            
  for (t in 1:M){
    log_y_obs[t] ~ normal(log_y_sim[t], sigma);
  }
  // priors
  theta ~ normal(0, 5);
  m0 ~ normal(300, 300);
}

generated quantities {
  real<lower=0> prod_theta2_m0;
  real<lower=0> prod_theta2_scale;
  real<lower=0> prod_m0_scale;
  real<lower=0> prod_theta2_m0_scale;
  vector[M] x2_sim;
  prod_theta2_m0 = theta[2] * m0;
  prod_theta2_scale = theta[2] * scale;
  prod_m0_scale = m0 * scale;
  prod_theta2_m0_scale = theta[2] * m0 * scale;
  x2_sim = ode_sol(t0, time_points, M, m0, theta);
}
