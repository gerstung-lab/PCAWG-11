data {
  int<lower=0> n;
  int<lower=0> p;          // number of data points
  vector<lower=0>[n] y;    // obs
  matrix[n,p] x;
  matrix[n,p] t;
}

parameters {
  real<lower=0> sigma; 
  real<lower=0> tau;
  real<lower=0> upsilon;
  real<lower=0> phi;
  real<lower=0> chi;
  vector<lower=0>[p] alpha;
  vector<lower=0>[p] beta;
  vector<lower=0>[p] delta;
  vector<lower=0>[n] mu; 
}

transformed parameters {
  vector[n] nu;
  vector[n] lambda;
  vector[n] epsilon;
  nu = x * (beta .* delta);
  lambda = mu + t * alpha;
  epsilon = t * delta;
}

model {
 beta ~ gamma(tau, upsilon);
 alpha ~ gamma(phi, chi);
 mu ~ gamma(nu, epsilon);
 y ~ normal(lambda, sigma); 
}
