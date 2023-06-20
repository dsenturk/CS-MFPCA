CS_MFPCA_simulation <- function(ngrid = 32, # Number of regions (scalar)
                                N = 100, # Number of individuals (scalar)
                                RetentPct = 0.5, # Percentage that you want to reserve in the final analysis (scalar between 0 and 1)
                                sigma = 0.1 # Measurement error variance (scalar)
){
#############################################################################
## Description: Function for simulating one data set under the simulation design described
##              in Section 6.
## Args: see above
## Definition: ngrid: number of total time points;  
##             M: number of dimensions in the original multivariate functional outcome X(t);
##             N: number of total individuals.
## Returns: list()
  ##     data.T: list of all true values, including:
  #           X0.T: underlying true original multivariate functional outcome in the reference dimension (matrix of dimension N*ngrid)
  #           X1.T: underlying true original multivariate functional outcome in the 1st dimension (matrix of dimension N*ngrid)
  #           X2.T: underlying true original multivariate functional outcome in the 2nd dimension (matrix of dimension N*ngrid) 
  #           X3.T: underlying true original multivariate functional outcome in the 3rd dimension (matrix of dimension N*ngrid)
  #           Z1.T: underlying true transformed multivariate functional outcome in the 1st dimension (matrix of dimension N*ngrid)
  #           Z2.T: underlying true transformed multivariate functional outcome in the 2nd dimension (matrix of dimension N*ngrid) 
  #           Z3.T: underlying true transformed multivariate functional outcome in the 3rd dimension (matrix of dimension N*ngrid)
  #           gridPoints: trial points (vector of length ngrid)
  #           psi1.True: underlying true multivariate eigenfunctions in the first dimension (matrix of dimension ngrid*3)
  #           psi2.True: underlying true multivariate eigenfunctions in the second dimension (matrix of dimension ngrid*3)
  #           psi3.True: underlying true multivariate eigenfunctions in the third dimension (matrix of dimension ngrid*3)
  #           v.True: underlying true multivariate eigen values (vector of length 3)
  ##     data.obs: list of observed values, including:
  #           X0: original multivariate functional outcome in the reference dimension (matrix of dimension N*ngrid)
  #           X1: original multivariate functional outcome in the 1st dimension (matrix of dimension N*ngrid)
  #           X2: original multivariate functional outcome in the 2nd dimension (matrix of dimension N*ngrid) 
  #           X3: original multivariate functional outcome in the 3rd dimension (matrix of dimension N*ngrid)
  #           Z1: transformed multivariate functional outcome in the 1st dimension (matrix of dimension N*ngrid)
  #           Z2: transformed multivariate functional outcome in the 2nd dimension (matrix of dimension N*ngrid) 
  #           Z3: transformed multivariate functional outcome in the 3rd dimension (matrix of dimension N*ngrid) 
  
################################################################################## 

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


# First, generate the multivariate functional outcome Z(t) with (M-1) dimensions
# Generate multivariate eigenfunctions in each dimension
fbs <- create.fourier.basis(rangeval = c(0, 3), nbasis = 13)
gridPoints <-  seq(0, 1, length.out = ngrid)

# Multivariate eigenfunction in the 1st dimension
phi1.True <- phi11 <- eval.basis(gridPoints, fbs)[,c(6, 7, 12)]
for(e in 1:3){ 
  coef <- sample(c(1,-1), 1, replace = TRUE)
  phi1.True[,e] <- coef * phi11[,e]
}

# Multivariate eigenfunction in the 2nd dimension 
x2 <- seq(1, 2, length.out = ngrid)
phi2.True <- phi22 <- eval.basis(x2, fbs)[,c(7, 12, 6)]
for(e in 1:3){ 
  coef <- sample(c(1,-1), 1, replace = TRUE)
  phi2.True[,e] <- coef * phi22[,e]
}

# Multivariate eigenfunction in the 3rd dimension 
x3 <- seq(2, 3, length.out = ngrid)
phi3.True <- phi33 <- eval.basis(x3, fbs)[,c(12, 6, 7)]
for(e in 1:3){ 
  coef <- sample(c(1,-1), 1, replace = TRUE)
  phi3.True[,e] <- coef * phi33[,e]
}

# Eigen scores in the multivariate functional outcome Z(t)
rho <- array(rep(NA, N*3), c(N, 3))
v.True <- rep(NA, 3)
for (i in 1:3){
  set.seed(6+i)
  v.True[i] <- (4 - i)/4
  rho[,i] <- rnorm(N, mean = 0, sd = sqrt(v.True[i]))
}


# Mean functions in three dimensions for Z(t)
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

# Measurement error variance
sigma <- sigma


# Define the location of removed trials for each individual
RetentNgrid <- ceiling(ngrid*RetentPct)
rm.cnt <- array(rep(NA, N*(ngrid-RetentNgrid)), c(N, (ngrid-RetentNgrid)))
for (i in 1:N){
  set.seed(1+i)
  rm.cnt[i,] <- sample(1:ngrid, (ngrid-RetentNgrid), replace = FALSE)
}



#######################################################
# Construct data set for three-dimensional multivariate 
# functional outcome Z(t)
#######################################################
# 1st dimension
Z1.0 <- Z1.0.T <- Z1 <- Z1.T <- matrix(0.0, nrow = N, ncol = ngrid)
eps1e <- matrix(rnorm(N * ngrid, 0, sqrt(sigma)), nrow = N)
for(i in 1:N){
  Z1[i,] <- Z1.0[i,] <- mu1(gridPoints) + rho[i,1] * phi1.True[,1] + 
                        rho[i,2] * phi1.True[,2] + rho[i,3] * phi1.True[,3] + 
                        eps1e[i,]
  Z1.T[i,] <- Z1.0.T[i,] <- mu1(gridPoints) + rho[i,1] * phi1.True[,1] + 
                            rho[i,2] * phi1.True[,2] + rho[i,3] * phi1.True[,3]
}

for (i in 1:N){
  Z1[i, c(rm.cnt[i,])] <- NA
  Z1.T[i, c(rm.cnt[i,])] <- NA
}

# 2nd dimension
Z2 <- Z2.T <- Z2.0 <- Z2.0.T <- matrix(0.0, nrow = N, ncol = ngrid)
eps2e <- matrix(rnorm(N * ngrid, 0, sqrt(sigma)), nrow = N)
for(i in 1:N){
  Z2[i,] <- Z2.0[i,] <-  mu2(gridPoints) + rho[i,1] * phi2.True[,1] + 
                         rho[i,2] * phi2.True[,2] + rho[i,3] * phi2.True[,3] +
                         eps2e[i,]
  Z2.T[i,] <- Z2.0.T[i,] <- mu2(gridPoints) + rho[i,1] * phi2.True[,1] + 
                            rho[i,2] * phi2.True[,2] + rho[i,3] * phi2.True[,3]
}

for (i in 1:N){
  Z2[i, c(rm.cnt[i,])] <- NA
  Z2.T[i, c(rm.cnt[i,])] <- NA
}

# 3rd dimension
Z3 <- Z3.T <- Z3.0 <- Z3.0.T <- matrix(0.0, nrow = N, ncol = ngrid)
eps3e <- matrix(rnorm(N * ngrid, 0, sqrt(sigma)), nrow = N)
for(i in 1:N){
  Z3[i,] <- Z3.0[i,] <- mu3(gridPoints) + rho[i,1] * phi3.True[,1] + 
                        rho[i,2] * phi3.True[,2] + rho[i,3] * phi3.True[,3] + 
                        eps3e[i,]
  Z3.T[i,] <- Z3.0.T[i,] <- mu3(gridPoints) + rho[i,1] * phi3.True[,1] + 
                            rho[i,2] * phi3.True[,2] + rho[i,3] * phi3.True[,3] 
}

for (i in 1:N){
  Z3[i, c(rm.cnt[i,])] <- NA
  Z3.T[i, c(rm.cnt[i,])] <- NA
}



###############################################################
# Construct data set for four-dimensional multivariate 
# functional outcome X(t)
###############################################################
# Underlying true multivariate functional outcome X(t)
X0.T <- X1.T <- X2.T <- X3.T <- Z1
for (i in 1:N){
  for (j in 1:(ngrid)){
    X0.T[i,j] <- 1/(1+exp(Z1.0.T[i,j])+exp(Z2.0.T[i,j]) + exp(Z3.0.T[i,j]))
    X1.T[i,j] <- exp(Z1.0.T[i,j])/(1+exp(Z1.0.T[i,j])+exp(Z2.0.T[i,j]) + exp(Z3.0.T[i,j]))
    X2.T[i,j] <- exp(Z2.0.T[i,j])/(1+exp(Z1.0.T[i,j])+exp(Z2.0.T[i,j]) + exp(Z3.0.T[i,j]))
    X3.T[i,j] <- exp(Z3.0.T[i,j])/(1+exp(Z1.0.T[i,j])+exp(Z2.0.T[i,j]) + exp(Z3.0.T[i,j]))
  }
}

# Observed multivariate functional outcome X(t) with noise
X0 <- X1 <- X2 <- X3 <- Z1
for (i in 1:N){
  for (j in 1:(ngrid)){
    X0[i,j] <- 1/(1+exp(Z1[i,j])+exp(Z2[i,j]) + exp(Z3[i,j]))
    X1[i,j] <- exp(Z1[i,j])/(1+exp(Z1[i,j])+exp(Z2[i,j]) + exp(Z3[i,j]))
    X2[i,j] <- exp(Z2[i,j])/(1+exp(Z1[i,j])+exp(Z2[i,j]) + exp(Z3[i,j]))
    X3[i,j] <- exp(Z3[i,j])/(1+exp(Z1[i,j])+exp(Z2[i,j]) + exp(Z3[i,j]))
  }
}

# Store all the true values and datasets
data.T <- list(X0.T = X0.T, X1.T = X1.T, X2.T = X2.T, X3.T = X3.T,
               Z1.T = Z1.T, Z2.T = Z2.T, Z3.T = Z3.T,
               gridPoints = gridPoints, phi1.True = phi1.True, 
               phi2.True = phi2.True, phi3.True = phi3.True,
               v.True = v.True)

data.obs <- list(X0 = X0, X1 = X1, X2 = X2, X3 = X3, 
                 Z1 = Z1, Z2 = Z2, Z3 = Z3)

# Generate the outputs
out <- list(data.T = data.T, data.obs = data.obs)

return(out)
}












