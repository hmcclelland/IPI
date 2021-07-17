#  -------------------------------------------------------
# Supplement to: "Statistical Uncertainty in Paleoclimate Proxy Reconstructions"
# Please cite:  McClelland et al. 2021., GRL.  DOI: 10.1029/2021GL092773
#  -------------------------------------------------------
#  This R script and accompanying files can be used to recreate Figures S2 A,B,D in McClelland et al. 2021. 
#  
#  correspondence to: Harry McClelland: hmcclelland@unimelb.edu.au (University of Melbourne, School of Earth Sciences)
#
#  To run this script, the following files need to be in the working directory: 
#   "Generate_data.R", "def_IPI_func", "mod_2stage_cal.stan", "mod_2stage_IP_1miltinewy.stan"

#  Recommended: Keep downloaded files together and set working directory to source file location.  



# Generate artificial data: 
source('Generate_data.R')

# Load artificial calibration data: 
Cal_data1 <- read.table("Cal_data.csv", header = TRUE, sep = ",", row.names = 1)
X_cal <- Cal_data1[,1]
Y_cal <- Cal_data1[,2]
New_y1 <- seq(min(Y_cal), max(Y_cal), length.out = 20)

# Load artificial paleo data: 
Series_data1 <- read.table("t_series_data.csv", header = TRUE, sep = ",", row.names = 1)
New_y2 <- Series_data1[,2]
ts <- Series_data1[,1]

# compile Bayesian models in stan 
library("rstan")
rstan_options(auto_write = TRUE)
mod_2step_cal <- stan_model('mod_2stage_cal.stan')
mod_2step_IP <- stan_model('mod_2stage_IP_1miltinewy.stan')

# Define IPI function: 
source('def_IPI_func.R')
# Inputs: Calibration X and Y values, and new Y values.  Confidence level is set tp 95% by default. 


#-------------------------------------------------------------------------------------------------
# Do Calculations  ( can substitute in own data here instead of using the artificial data set ) :
#-------------------------------------------------------------------------------------------------

# Calibration: 
Cal_IPI <- IPIfunc(X_cal, Y_cal, New_y1)

# Reconstruction: 
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
  
  Cqts   <- c((1-Clev/100)/2, 1 - (1-Clev/100)/2)
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
  
  
  
  
  