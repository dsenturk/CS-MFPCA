## CS_MFPCA_tutorial.R
#############################################################################
## Description: A step-by-step implementation of CS-MFPCA and the associated  
## procedures described in "Constrained multivariate functional principal  
## components analysis for novel outcomes in eye-tracking experiments".  
#############################################################################
## Functions implemented: 
## CS_MFPCA_simulation.R (Simulate a four-dimensional multivariate functional dataset), 
## CS_MFPCA_algorithm.R (Step 1-2 of the CS-MFPCA algorithm),
## CS_MFPCA_prediction.R (Step 3 of the CS-MFPCA algorithm).
#############################################################################
## Tutorial Outline:
## 1. Simulate four-dimensional multivariate functional outcome (CS_MFPCA_simulation.R)
## 2. Perform CS-MFPCA algorithm (CS_MFPCA_algorithm.R, CS_MFPCA_prediction.R)
## 3. Measurement of CS-MFPCA results
#############################################################################

# Install missing packages
list.of.packages <- c("refund", "fda", "funData", "refund", "caTools", "locpol", 
                      "htmltools", "fdapace", "abind")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 

# Load packages
library(funData)
library(MFPCA)
library(fda)
library(refund)
library(caTools)
library(locpol)
library(htmltools)
library(fdapace)
library(abind)

#############################################################################
# 1. Simulate the original four-dimensional multivariate functional outcome X(t)
#   as well as transformed three-dimensional multivariate functional outcome Z(t)
#############################################################################

# Simulate one dataset from the simulation design described 
data.G <- CS_MFPCA_simulation(ngrid = 16, N = 100, RetentPct = 1, sigma = 0.1)  
# CS_MFPCA_simulation.R (ngrid = 32 if you want in total 32 time points, 
                      # N = 100 if you want in total 100 individuals,
                      # RetentPct = 0.5 (RetentPct is the percentage that you want to reserve in the final analysis 
                                         # from the original generated time points) if you want to reserve 50% of 
                                         # the `ngrid` time points.
                                         # Note that when RetentPct = 1 is the regular design and 
                                         # RetentPct = 0.5 is the irregular design in Section 6.
                      # sigma = 0.1 where 0.1 can be any scalar representing the setup error variance)

# Data array used for estimation
data.obs <- data.G$data.obs

# True underlying data list 
data.T <- data.G$data.T

#############################################################################
# 2. Perform CS-MFPCA algorithm
#############################################################################

# NOTE: Performing CS-MFPCA algorithm Steps 2 

MFPCAout <- CS_MFPCA_algorithm(data.obs)  # CS_MFPCA_algorithm.R



# NOTE: Performing CS-MFPCA estimation steps 3

Predout <- CS_MFPCA_prediction(MFPCAout) # CS_MFPCA_prediction.R  


#############################################################################
# 3. Measurement of CS-MFPCA results
#############################################################################  
# Define the total number of time points, specific time points and total individuals
gridPoints <- data.T$gridPoints
ngrid <- length(gridPoints)
N <- nrow(data.T$X0.T)

# True mean functions for Z(t)
# 1st dimension
mu1 <- function(x){
  return(0.4*x + 2)
}

# 2nd dimension
mu2 <- function(x){
  return(0.8*x + 2)
}

# 3rd dimension
mu3 <- function(x){
  return(0.2*x + 2)
}


## Define RSE function to calculate the prediction error 
RSE.fun <- function(x,y){
  err1 <- trapz(gridPoints, (x-y)^2)
  err2 <- trapz(gridPoints, (x+y)^2)
  return(c(min(err1, err2)))
}

# RSE for predictions
RSE.X0.pre <- rep(0, N)
RSE.X1.pre <- rep(0, N)
RSE.X2.pre <- rep(0, N)
RSE.X3.pre <- rep(0, N)
for (i in 1:N){
  RSE.X0.pre[i] <- RSE.fun(data.T$X0.T[i,], Predout$X0.est[i,])
  RSE.X1.pre[i] <- RSE.fun(data.T$X1.T[i,], Predout$X1.est[i,])
  RSE.X2.pre[i] <- RSE.fun(data.T$X2.T[i,], Predout$X2.est[i,])
  RSE.X3.pre[i] <- RSE.fun(data.T$X3.T[i,], Predout$X3.est[i,])
}
RSE.X0 <- mean(RSE.X0.pre)
RSE.X1 <- mean(RSE.X1.pre)
RSE.X2 <- mean(RSE.X2.pre)
RSE.X3 <- mean(RSE.X3.pre)

# RSE for multivariate eigenfunctions
RSE.phi1 <- RSE.phi2 <- RSE.phi3 <- rep(NA, 3)
for (i in 1:3){
  RSE.phi1[i] <- RSE.fun(data.T$phi1.True[,i], MFPCAout$phi1.est[,i])
  RSE.phi2[i] <- RSE.fun(data.T$phi2.True[,i], MFPCAout$phi2.est[,i])
  RSE.phi3[i] <- RSE.fun(data.T$phi3.True[,i], MFPCAout$phi3.est[,i])
}

# RSE for mean functions
RSE.mu <- rep(NA, 3)
for (i in 1:3){
  RSE.mu[i] <- RSE.fun(mu1(gridPoints), MFPCAout$mu1.est)
  RSE.mu[i] <- RSE.fun(mu2(gridPoints), MFPCAout$mu2.est)
  RSE.mu[i] <- RSE.fun(mu3(gridPoints), MFPCAout$mu3.est)
}

# MSE for eigenvalues
(data.T$v.True - MFPCAout$v.est)^2/data.T$v.True^2



##############################################################
# Plot estimates of original multivariate trajectory predictions 
# in four dimensions
##############################################################

# Predicted trajectory for a specific individual
# Subject ID for plots
ID <- 1
# True trajectory for the `ID`th individual
traj0.T <- data.T$X0.T[ID,] # reference dimension
traj1.T <- data.T$X1.T[ID,] # 1st dimension
traj2.T <- data.T$X2.T[ID,] # 2nd dimension
traj3.T <- data.T$X3.T[ID,] # 3rd dimension
# Predicted trajectory for the R.ID`th individual
traj0.Est <- Predout$X0.est[ID,] # reference dimension
traj1.Est <- Predout$X1.est[ID,] # 1st dimension
traj2.Est <- Predout$X2.est[ID,] # 1st dimension
traj3.Est <- Predout$X3.est[ID,] # 1st dimension

par(mfrow = c(2,2))
# reference dimension
ylim1 <- min(traj0.T, traj0.Est)
ylim2 <- max(traj0.T, traj0.Est)
plot(gridPoints, traj0.T, "l", lwd = 2, ylim = c(ylim1,ylim2), xaxs = "i", 
     main = "Reference dimension", cex.main = 2, xlab = "", ylab = "", col = "black")
title(xlab = "t", ylab = "Trajectory", line = 2, cex.lab = 1.6)
lines(gridPoints, traj0.Est, col = "grey52", lwd = 2, lty = 1)
legend("topright", legend = c("True", "Predicted"), col = c("black", "grey52"), 
       lty = c(1,1), cex = 0.85)

# 1st dimension
ylim1 <- min(traj1.T, traj1.Est)
ylim2 <- max(traj1.T, traj1.Est)
plot(gridPoints, traj1.T, "l", lwd = 2, ylim = c(ylim1,ylim2), xaxs = "i",
     main = "1st dimension", cex.main = 2, xlab = "", ylab = "", col = "black")
title(xlab = "t",ylab = "Trajectory", line = 2, cex.lab = 1.6)
lines(gridPoints, traj1.Est, col = "grey52", lwd = 2, lty = 1)

# 2nd dimension
ylim1 <- min(traj2.T, traj2.Est)
ylim2 <- max(traj2.T, traj2.Est)
plot(gridPoints, traj2.T, "l", lwd = 2, ylim = c(ylim1,ylim2), xaxs = "i", 
     main = "2nd dimension", cex.main = 2, xlab = "", ylab = "", col = "black")
title(xlab = "t", ylab = "Trajectory", line = 2, cex.lab = 1.6)
lines(gridPoints, traj2.Est, col = "grey52", lwd = 2, lty = 1)


# 3rd dimension
ylim1 <- min(traj3.T, traj3.Est)
ylim2 <- max(traj3.T, traj3.Est)
plot(gridPoints, traj3.T, "l", lwd = 2, ylim = c(ylim1,ylim2), xaxs = "i",
     main = "3rd dimension", cex.main = 2, xlab = "", ylab = "", col = "black")
title(xlab = "t",ylab = "Trajectory", line = 2, cex.lab = 1.6)
lines(gridPoints, traj3.Est, col = "grey52", lwd = 2, lty = 1)

