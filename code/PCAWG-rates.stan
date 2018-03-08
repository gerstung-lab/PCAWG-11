data {
  int<lower=0> n;		// numbr of observations
  int<lower=0> p;          // number of types
  vector<lower=0>[n] y;    // mutations
  matrix[n,p] x;		// age
  matrix[n,p] t;		// tumour type
}

parameters {
  real<lower=0> sigma; // const. variance 
  real<lower=0> tau; // time-dep variance
  real<lower=0> alpha;      // alpha for slope
  real<lower=0> beta;  // beta for slope
  real<lower=0> gamma;
  real<lower=0> delta;
  vector<lower=0>[p] offset;
  vector<lower=0>[p] slope;
}

transformed parameters {
  vector[n] mu;
  vector[n] nu;
  vector[p] ones;
  ones = rep_vector(1, p);
  mu = x * slope + t * offset;
  nu = sqrt( square(x * ones)  * tau^2 + sigma^2);
}

model {
 slope ~ gamma(alpha, beta);
 offset ~ gamma(gamma, delta);
 y ~ normal(mu, nu); 
}
