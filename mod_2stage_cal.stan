// Bayesian Inverse Prediction code for Stan.  
// Use with accommpanying R code IPI_Bayesian.R
// Harry McClelland -- 3/4/21
data {
  int<lower=0> n;       // no samples
  vector[n] x;          // calibration data of length N
  vector[n] y;  
}

parameters {
  real B0;              // No prior specified so assumes uniform prior. 
  real B1;
  real<lower=0> sig;    // Positive uniform prior. 
}

model {
  y ~ normal(B0 + B1 * x, sig);
}
