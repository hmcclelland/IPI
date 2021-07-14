// Bayesian Inverse Prediction code for Stan.  
// Use with accommpanying R code IPI_Bayesian.R
// Harry McClelland -- 3/4/21
data {
  int<lower=0> n;       // no samples
  vector[n] x;          // calibration data of length N
  vector[n] y;  
  
  real pr_B0_mu;
  real<lower=0> pr_B0_s;
  real pr_B1_mu;
  real<lower=0> pr_B1_s;
  real<lower=0> pr_sig;
}

parameters {
  real B0;              // No prior specified so assumes uniform prior. 
  real B1;
  real<lower=0> sig;    // Positive uniform prior. 
}

model {
  y ~ normal(B0 + B1 * x, sig);
  B0 ~ normal(pr_B0_mu, pr_B0_s);
  B1 ~ normal(pr_B1_mu, pr_B1_s);
  sig ~ exponential(pr_sig);
}
