data {
  int<lower=0> n;
  int<lower=0> p;          // number of data points
  vector[n] y;             // obs
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
}

transformed parameters {
  vector[n] mu;
  mu = x * beta +  t * alpha;
}

model {
 beta ~ gamma(tau, upsilon);
 alpha ~ gamma(phi, chi);
 y ~ normal(mu, sigma); 
}
