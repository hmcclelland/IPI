# Define function for calculating IPI: 

#  The function IPIfunc requires X and Y calibration values (vectors of equal length, n) and a vector of new Y0 values (length, m), 
#  with an optional specification of confidence level in % (defaults to 95%)
#
#  Usage:  IPIfunc(X_cal (1xn vector), Y_cal (1xn vector), Y_new (1xm vector), Clev (defaults to 95%))
#
#  The output is a table of estimated X, plus upper and lower values for IPIs calculated with 4 different appraoches. 
#
#  Output: mx10 array with column headings: 
#          "Y_new", "X_new_hat", 
#           "Simple_lower_X", "Simple_upper_X", 
#           "Nonpara_lower_X", "Nonpara_upper_X", 
#           "Fiducial_lower_X", "Fiducial_upper_X", 
#           "Bayesian_lower_X", "Bayesian_upper_X")
#
#  -------------------------------------------------------

IPIfunc <- function(X_cal, Y_cal, Y_new, Clev = 95){
  set.seed(26442)
  
  # load calibration data set:
  Xs     <- X_cal          #  X values are assumed to be known exactly. 
  Ys     <- Y_cal                   #  Y is a stochastic variable. 
  nCal   <- length(Xs)
  Cqts   <- c((1-Clev/100)/2, 1 - (1-Clev/100)/2)      #  Confidence quantiles (ie 95% -> 0.025 & 0.975)
  
  
  #-------------------------------------------------------------------------------------------------
  #-------     Fit OLS linear model to data set     ------------------------------------------------ 
  #-------------------------------------------------------------------------------------------------
  
  lm1        <- lm(Ys ~ Xs)                 # fit model to calibration data 
  lm1_Bs     <- lm1$coefficients             
  B0_hat     <- lm1_Bs[1]
  B1_hat     <- lm1_Bs[2]
  X_new_hat  <- (Y_new - B0_hat) / B1_hat
  
  
  #-------------------------------------------------------------------------------------------------
  #--------     IPI: rotated sd residuals approach  ("Simple method" in paper)  --------------------
  #-------------------------------------------------------------------------------------------------
  
  Tval   <- qt(Cqts, df=(nCal-2))       # calculate t-values at given confidence level
  Yhat   <- B0_hat + Xs*B1_hat 
  ysd    <- sd(Ys - Yhat)
  Xe1    <- (ysd * Tval) / B1_hat
  Xe1L   <- X_new_hat + Xe1[1]
  Xe1U   <- X_new_hat + Xe1[2]
  
  
  #-------------------------------------------------------------------------------------------------
  #--------     IPI: non-parametric appraoch      (suggested Method B in paper)  --------------------
  #-------------------------------------------------------------------------------------------------
  
  Xhat   <- (Ys - B0_hat) / B1_hat 
  Xres   <- Xs - Xhat
  Xe2    <- quantile(Xres, Cqts)        # 2.5th and 97.5th percentiles of residuals in X axis. 
  Xe2L   <- X_new_hat + Xe2[1]
  Xe2U   <- X_new_hat + Xe2[2]
  
  
  #-------------------------------------------------------------------------------------------------
  #-------      IPI: 'Fiducial' approach:    From Draper & Smith 1998 (3rd edition) Eq. 3.2.6 (p.84)------
  #-------------------------------------------------------------------------------------------------
  
  tv1   <- Tval[2]  # +ve T value
  SSyy  <-  sum((Ys - Yhat)^2)
  SSxx  <-  sum((Xs - Xhat)^2)
  s2    <-  SSyy / (nCal-2)
  Xbar  <- mean(Xs)
  Ybar  <- mean(Ys)
  q     <- 1                                  # q is number of repeats. if actual Y is known exactly, set Q as inf. (not relevant for proxies: "what is a true replicate")
  
  Qe_A <- B1_hat^2 - tv1^2 * s2 / SSxx
  Qe_B <- 2*(Xbar * tv1^2 * s2 / SSxx - B1_hat^2 * X_new_hat)
  Qe_C <- B1_hat^2 * X_new_hat^2 - tv1^2 * s2 *(1/q + 1/nCal) - tv1^2 * s2 * Xbar^2 / SSxx
  
  Xe3U <- ((-Qe_B) + (Qe_B^2 - 4*Qe_A*Qe_C)^(1/2))/ (2*Qe_A)
  Xe3L <- ((-Qe_B) - (Qe_B^2 - 4*Qe_A*Qe_C)^(1/2))/ (2*Qe_A)
  
  
  #-------------------------------------------------------------------------------------------------
  #---------------  IPI: Bayesian     -------------------------
  #-------------------------------------------------------------------------------------------------
  
  # Iterations of MCMC for each step. 
  iterations_2S_cal <- 2000
  iterations_2S_IP  <- 1000
  
  # NB: Priors for all regression parameters and for new x vaklues are flat (sigma is strictly +ve) 
  
  #----------  2 Step model:  Calibration  -------------------------------------
  
  # pass data to stan and run model
  options(mc.cores = parallel::detectCores())
  fit1 <- sampling(mod_2step_cal, list(n=nCal, x=Xs, y=Ys), 
                   iter = iterations_2S_cal, warmup = floor(iterations_2S_cal / 2), chains = 4, control = list(adapt_delta = 0.90, max_treedepth = 10))
  # diagnose
  print(fit1)
  params1 <- extract(fit1)
  
  fitted_B0 <- params1$B0
  fitted_B1 <- params1$B1
  fitted_sig <-params1$sig
  n_posts1 <- length(fitted_B0)
  
  fit_summary <- summary(fit1)
  summarytable <- fit_summary$summary
  Rhats <- summarytable[,"Rhat"]
  calibOK <- Rhats[1] & Rhats[2] & Rhats[3] < 1.1   # Make sure calibration solutions converged
  
  Nnew   <- length(Y_new)
  Xe4 <- matrix(nrow = Nnew, ncol = 2)
  
  if(calibOK){
    
    #----------  2 Step model:  Inverse Prediction  ---------------------------
    
    # pass data to stan and run model
    options(mc.cores = parallel::detectCores())
    fit2 <- sampling(mod_2step_IP, list(n_new=Nnew, y_new = Y_new, posts = n_posts1,
                                        B0 = fitted_B0, B1 = fitted_B1, sig = fitted_sig), 
                     iter = iterations_2S_IP, warmup = floor(iterations_2S_IP / 2), chains = 4, control = list(adapt_delta = 0.90, max_treedepth = 10))
    # diagnose
    print(fit2)
    
    params2 <- extract(fit2)
    Xouts2 <- params2$x_new
    Xdims2 <-dim(Xouts2)
    for(xi in 1: Xdims2[2]){
      xis <- as.vector(Xouts2[,xi,])
      Xe4[xi,] <- quantile(xis, Cqts)
    }
  }
  
  Xe4L <- Xe4[,1]
  Xe4U <- Xe4[,2]
  
  
  
  #-------------------------------------------------------------------------------------------------
  #---------------  Output all     -------------------------
  #-------------------------------------------------------------------------------------------------
  
  output <- cbind(Y_new, X_new_hat, Xe1L, Xe1U, Xe2L, Xe2U, Xe3L, Xe3U, Xe4L, Xe4U)
  colnames(output) <- c("Y_new", "X_new_hat", "Simple_lower_X", "Simple_upper_X", "Nonpara_lower_X", "Nonpara_upper_X", 
                        "Fiducial_lower_X", "Fiducial_upper_X", "Bayesian_lower_X", "Bayesian_upper_X")
  return(output)
  
  #-------------------------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------------
}

