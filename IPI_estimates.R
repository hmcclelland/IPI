#  -------------------------------------------------------
# Supplement to: "Statistical Uncertainty in Paleoclimate Proxy Reconstructions"
# Please cite:  McClelland et al. 2021., GRL.  DOI: 10.1029/2021GL092773
#  -------------------------------------------------------
#  This R script and accompanying files can be used to recreate Figures S2 A,B,D in McClelland et al. 2021. 
#  
#  correspondence to: Harry McClelland: hmcclelland@unimelb.edu.au (University of Melbourne, School of Earth Sciences)
#
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


# Define function for calculating IPI: 
# compile Bayesian model in stan 
library("rstan")
rstan_options(auto_write = TRUE)
mod_2step_cal <- stan_model('mod_2stage_cal.stan')
mod_2step_IP <- stan_model('mod_2stage_IP_1miltinewy.stan')

#
# Inputs: Calibration X and Y values, and new Y values.  Confidence level is set tp 95% by default. 
#
IPIfunc <- function(X_cal, Y_cal, Y_new, Clev = 95){
  #set.seed(26441)
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
  
  # calibration model priors:
  pr_B0_mu = mean(Ys)
  pr_B0_s = 2*sd(Ys)    # NB:  need to check that priors are not strongly informative and reasonable for the data
  pr_B1_mu = 0          # 0
  pr_B1_s = 10
  pr_sig = 0.2
  
  # priors for new x:
  px_new_hat <- 0
  px_new_err <- 5
  
  # Iterations of MCMC for each step. 
  iters_2S_cal <- 500
  iters_2S_IP  <- 100
  
  #----------  2 Step model:  Calibration  -------------------------------------
  
  # pass data to stan and run model
  options(mc.cores = parallel::detectCores())
  fit1 <- sampling(mod_2step_cal, list(n=nCal, x=Xs, y=Ys, 
                                       pr_B0_mu = pr_B0_mu, pr_B0_s = pr_B0_s,
                                       pr_B1_mu = pr_B1_mu, pr_B1_s = pr_B1_s,
                                       pr_sig = pr_sig), 
                   iter = iters_2S_cal, chains = 4, control = list(adapt_delta = 0.90, max_treedepth = 10))
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
                                        B0 = fitted_B0, B1 = fitted_B1, sig = fitted_sig,
                                        px_new_hat = rep(px_new_hat, n_posts1), px_new_err = rep(px_new_err, n_posts1)), 
                     iter = iters_2S_IP, chains = 4, control = list(adapt_delta = 0.90, max_treedepth = 10))
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
  








#-------------------------------------------------------------------------------------------------
# Do Calculations:
#-------------------------------------------------------------------------------------------------


# Calibration 
Cal_data1 <- read.table("Cal_data.csv", header = TRUE, sep = ",", row.names = 1)
X_cal <- Cal_data1[,1]
Y_cal <- Cal_data1[,2]
New_y1 <- seq(min(Y_cal), max(Y_cal), length.out = 20)
Cal_IPI <- IPIfunc(X_cal, Y_cal, New_y1)

# Reconstruction
Series_data1 <- read.table("t_series_data.csv", header = TRUE, sep = ",", row.names = 1)
New_y2 <- Series_data1[,2]
ts <- Series_data1[,1]
Series_IPI <- IPIfunc(X_cal, Y_cal, New_y2)


#-------------------------------------------------------------------------------------------------
# Plot:
#-------------------------------------------------------------------------------------------------

pdf("Example_IPI_plots_A.pdf", width = 6, height = 6)
par(mfrow = c(1,1))

COLS        <- c("black", "#009E73", "#0072B2","#D55E00",  "#CC79A7", "black")  # ---->  plotting: Cols (1:5)
thicknesses <- c(4, 1.5, 1.5, 1.5, 1.5, 2)   # -------------------------------------->  plotting: line thicknesses (1:5)
linetypes   <- c(1, 2, 3, 4, 1, 1)   # --------------------------------------> plotting: line types (1:5)



# DATA: 
plot(X_cal, Y_cal, ylim=c(min(Cal_IPI[,1]), max(Cal_IPI[,1])), 
     xlim=c(min(Cal_IPI[,2:10]), max(Cal_IPI[,2:10])), pch = 21, bg = 'grey', col='white', 
       xlab = "Environmental variable, E", ylab = "Proxy variable, P")
# Regression line:
P = 1                                 # Plot number 
lines(Cal_IPI[,2], Cal_IPI[,1], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])

# IPI ESTIMATES: 
  # plot:
  P = 2                                 # Plot number 
  lines(Cal_IPI[,3], Cal_IPI[,1], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  lines(Cal_IPI[,4], Cal_IPI[,1], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  # plot:
  P = 3                                 # Plot number 
  lines(Cal_IPI[,5], Cal_IPI[,1], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  lines(Cal_IPI[,6], Cal_IPI[,1], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  # plot:
  P = 4                                 # Plot number 
  lines(Cal_IPI[,7], Cal_IPI[,1], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  lines(Cal_IPI[,8], Cal_IPI[,1], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  #plot:
  P = 5                                 # Plot number 
  lines(Cal_IPI[,9], Cal_IPI[,1], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  lines(Cal_IPI[,10], Cal_IPI[,1], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  
  
  legend("topleft", 
         legend = c("Calibration data", "OLS fitted line", "95% IPI: Simple", "95% IPI: Non-parametric", 
                    "95% IPI: Fiducial", "95% IPI: Bayesian"), pch = c(21, NA,NA,NA,NA, NA), pt.bg = c("grey", NA,NA,NA,NA, NA),
         col = c("white", COLS), lwd = c(NA, thicknesses), lty = c(NA, linetypes), cex=0.8, bty="n")
  
  CDF_AT<- 0.55
  abline(h = CDF_AT)
  text(0, 0.57, labels = "CDF (Fig D)", cex=0.9)
 
  axislims<-par("usr")
  mtext("A", line=1, at = axislims[1], cex=1.1)
  
  dev.off()
  
  # Reconstruction 
  

  
  
  # Plotting -- 
  pdf("Example_IPI_plots_B.pdf", width = 6, height = 6)
  par(mfrow = c(1,1))
  
  COLS        <- c("black", "#009E73", "#0072B2","#D55E00",  "#CC79A7")  # ---->  plotting: Cols (1:5)
  thicknesses <- c(4, 1.5, 1.5, 1.5, 1.5)   # -------------------------------------->  plotting: line thicknesses (1:5)
  linetypes   <- c(1, 2, 3, 4, 1)   # --------------------------------------> plotting: line types (1:5)
  
  plot(ts, Series_IPI[,2], ylim=c(min(Series_IPI[,2:10]), max(Series_IPI[,2:10])), pch = 21, bg = 'grey', col='white', cex=0.1,
       xlab = "time", ylab = "Reconstructed E")
  # plot:
  P = 1                                 # Plot number 
  lines(ts, Series_IPI[,2], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  # plot:
  P = 2                                 # Plot number 
  lines(ts, Series_IPI[,3], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  lines(ts, Series_IPI[,4], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  # plot:
  P = 3                                 # Plot number 
  lines(ts, Series_IPI[,5], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  lines(ts, Series_IPI[,6], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  # plot:
  P = 4                                 # Plot number 
  lines(ts, Series_IPI[,7], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  lines(ts, Series_IPI[,8], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  #plot:
  P = 5                                 # Plot number 
  lines(ts, Series_IPI[,9], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  lines(ts, Series_IPI[,10], col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  
  
  legend("topright", 
         legend = c("Estimate", "95% IPI: Simple", "95% IPI: Non-parametric", 
                    "95% IPI: Fiducial", "95% IPI: Bayesian"), pch = c(NA,NA,NA,NA, NA), pt.bg = c(NA,NA,NA,NA, NA),
         col = c(COLS), lwd = c(thicknesses), lty = c(linetypes), cex=0.8, bty="n")
  
  axislims<-par("usr")
  mtext("B", line=1, at = axislims[1], cex=1.1)
  
  dev.off()
  
  
  
  
  #-------------------------------------------------------------------------------------------------
  # MC Simulations:
  #-------------------------------------------------------------------------------------------------
  
  pdf("Example_IPI_plots_D.pdf", width = 6, height = 6)
  
  # MC simulation CDF plot
  par(mfrow = c(1,1))
  
  MC_data <- read.table("MCout.csv", header = TRUE, sep = ",", row.names = 1) 
  sortedX <- sort(MC_data$x)
  e_cdf <- 1:length(sortedX) / length(sortedX)
  plot(sortedX, e_cdf, type = "s", xlab = expression(paste(Delta, "E")), ylab = "Monte Carlo CDF")
  Lecd<- sortedX[which(e_cdf >= Cqts[1])[1]]
  Uecd<- sortedX[which(e_cdf >= Cqts[2])[1]]
  # "Offset from X"[0]
  polygon(c(sortedX, max(sortedX), min(sortedX)), c(e_cdf, 0, 0), col=rgb(0.9, 0.9, 0.9), border = NA)
  abline(h=c(0.025, 0.975))
  
  Xappx <- approx(Cal_IPI[,1], Cal_IPI[,2], CDF_AT)$y
  XappxU <- approx(Cal_IPI[,1], Cal_IPI[,4], CDF_AT)$y
  XappxL <- approx(Cal_IPI[,1], Cal_IPI[,3], CDF_AT)$y
  XU <- XappxU - Xappx
  XL <- XappxL - Xappx
  P = 2                                 # Plot number 
  abline(v=c(XU, XL), col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  XappxU <- approx(Cal_IPI[,1], Cal_IPI[,6], CDF_AT)$y
  XappxL <- approx(Cal_IPI[,1], Cal_IPI[,5], CDF_AT)$y
  XU <- XappxU - Xappx
  XL <- XappxL - Xappx
  P = 3                                 # Plot number 
  abline(v=c(XU, XL), col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  XappxU <- approx(Cal_IPI[,1], Cal_IPI[,8], CDF_AT)$y
  XappxL <- approx(Cal_IPI[,1], Cal_IPI[,7], CDF_AT)$y
  XU <- XappxU - Xappx
  XL <- XappxL - Xappx
  P = 4                                 # Plot number 
  abline(v=c(XU, XL), col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  XappxU <- approx(Cal_IPI[,1], Cal_IPI[,10], CDF_AT)$y
  XappxL <- approx(Cal_IPI[,1], Cal_IPI[,9], CDF_AT)$y
  XU <- XappxU - Xappx
  XL <- XappxL - Xappx
  P = 5                                 # Plot number 
  abline(v=c(XU, XL), col = COLS[P],  lty = linetypes[P], lwd=thicknesses[P])
  
  legend("left", 
         # legend = c("95% IPI: Simple", "95% IPI: Non-parametric", 
         #            "95% IPI: Fiducial", "95% IPI: Bayesian"),
         legend = c("Simple", "Non-parametric", 
                    "Fiducial", "Bayesian"),
         col = c(COLS[2:5]), lwd = c(thicknesses[2:5]), lty = c(linetypes[2:5]), cex=0.8, bty="n")
  text(0, Cqts[1]+0.04, expression(paste("2.5"^"th", "percentile")), cex=0.8, col = 1)
  text(0, Cqts[2]-0.04, expression(paste("97.5"^"th", "percentile")), cex=0.9, col = 1)
  
  axislims<-par("usr")
  
  mtext("D", line=1, at = axislims[1], cex=1.1)
  
  dev.off()
  
  
  
  
  