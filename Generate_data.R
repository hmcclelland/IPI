#  -------------------------------------------------------
# Supplement to McClelland et al. 20##...
# Please cite the following when using:  
#  -------------------------------------------------------
#  This R script and accompanying files can be used to recreate Figures S2 A,B,D in McClelland et al. 20##. 
#  See accompanying readme file for details of use. 
#  correspondence to: Harry McClelland: hmcclelland@unimelb.edu.au (University of Melbourne, School of Earth Sciences)
#  -------------------------------------------------------


# Generate synthetic calibration data: 
set.seed(26473)


Beta0  <- 0.1       #-----------------------------------------------------> 'True' intercept 
Beta1  <- 0.8       #-----------------------------------------------------> 'True' slope
sigma  <- 0.07      #-----------------------------------------------------> 'True' sd of noise in Y axis
nCal   <- 100        #-----------------------------------------------------> n points in calibration data set. 
Xrange <- c(0,1)    #-----------------------------------------------------> min and max X for calibration. 

Xs     <- runif(nCal, Xrange[1], Xrange[2])          #  X values are assumed to be known exactly. 
noise  <- rnorm(nCal, mean = 0, sd = sigma)          #  noise is assumed normal
Ys     <- Beta0 + Xs*Beta1 + noise                   #  Y is a stochastic variable. 
Cal_data <- cbind(Xs, Ys)
write.table(Cal_data, "Cal_data.csv", sep = ",", qmethod = "double")


# Generate synthetic reconstruction data: 

xcal <- 1 # noise in time series relative to noise in calibration 

t     <- seq(1, 30)   # Time points
t_series_noise  <- rnorm(length(t), mean = 0, sd = xcal*sigma)          #  noise is assumed normal
t_series_Xs     <- 3*sigma*sin(t/5)  + mean(Xs)
t_series_Ys     <- Beta0 + t_series_Xs*Beta1 + t_series_noise        #  Y is a stochastic variable. 
t_series_data <- cbind(t, t_series_Ys)
write.table(t_series_data, "t_series_data.csv", sep = ",", qmethod = "double")






#
# The following code uses the "true" parameter values to simulate 
#
# Test with Monte Carlo simulation?  Run_MC  -->  1 = y,  0 = n?

Run_MC = 1

if(Run_MC == 1){
  
  Cal_data1 <- read.table("Cal_data.csv", header = TRUE, sep = ",", row.names = 1)
  X_cal <- Cal_data1[,1]
  Y_cal <- Cal_data1[,2]
  Clev <- 95
  Cqts   <- c((1-Clev/100)/2, 1 - (1-Clev/100)/2)      #  Confidence quantiles (ie 95% -> 0.025 & 0.975)
  
  TEST_CDF_AT <- 0.55
  CDF_Window <- 0.05
    
  nMC <- 100000 #  no. MC runs
  Nnew <- 300
  
  Xnews  <- seq(-3, 4, length.out=Nnew)             # "true"  new values of X

  lm1        <- lm(Y_cal ~ X_cal)                 # fit model to calibration data 
  lm1_Bs     <- lm1$coefficients             
  B0_hat     <- lm1_Bs[1]                         # OLS estimate of intercept and slope based on Cal_data.csv
  B1_hat     <- lm1_Bs[2]
  
  Y_out  <- matrix(nrow = Nnew, ncol = nMC)          #  preallocated vector for binned MC Y
  X_out  <- matrix(nrow = Nnew, ncol = nMC)         #   preallocated vector for binned MC Error in Predicted X 
  
  Xnews  <- seq(-3, 4, length.out=Nnew) 
  
  for(m in seq(1,nMC)){                                                 #  MC simulation 
    Ynew <-  Xnews*Beta1 + Beta0 + rnorm(Nnew, mean = 0, sd = sigma)    #  Generate stochastic new Ys from 'True' new Xs 
    Xest <-  (Ynew - B0_hat) / B1_hat                                         #  Estimate Xs based on new Ys using calibration estimates of B0 and B1
    Xerr <-  Xest - Xnews                                               #  Error in estimate of X compared with 'true' X
    Y_out[,m] <- Ynew                                        
    X_out[,m] <- Xerr
  }
  
  Yout = c(Y_out);            # Create vector from matrix
  Xout = c(X_out);            # Create vector from matrix

  index <- (Yout > (TEST_CDF_AT - CDF_Window/2) & Yout < (TEST_CDF_AT + CDF_Window/2));        # find index of those values in bin, b
  MCout <- Xout[index];

  write.table(MCout, "MCout.csv", sep=",")
}
 

 