// Bayesian Inverse Prediction code for Stan.  
// Use with accommpanying R code IPI_Bayesian.R
// Harry McClelland -- 3/4/21
data {
  int<lower=0> n_new;   // No 'new' y measurements to infer x_new from
  vector[n_new] y_new;

  int<lower=0> posts;
  vector[posts] B0;              
  vector[posts] B1;
  vector[posts] sig;   
}

parameters {
  matrix[n_new, posts] x_new;   # NB: flat priors
}

model {
  vector[posts] y_new_hat; 
  for(i in 1:n_new){
    y_new_hat = B0 + B1 .* x_new[i,]';
    y_new[i] ~ normal(y_new_hat, sig);
}
}
