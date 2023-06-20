CS_MFPCA_prediction <- function(MFPCAout # Output from function CS_MFPCA_MFPCA including the following: 
                               ##              Z1.est: estimates of transformed multivariate functional outcome in the 1st dimension (matrix of dimension N*ngrid)
                               ##              Z2.est: estimates of transformed multivariate functional outcome in the 2nd dimension (matrix of dimension N*ngrid)
                               ##              Z3.est: estimates of transformed multivariate functional outcome in the 3rd dimension (matrix of dimension N*ngrid)
                               ##              phi1.est: estimated multivariate eigenfunctions in the 1st dimension (matrix of dimension ngrid*3)
                               ##              phi2.est: estimated multivariate eigenfunctions in the 2nd dimension (matrix of dimension ngrid*3)
                               ##              phi3.est: estimated multivariate eigenfunctions in the 3rd dimension (matrix of dimension ngrid*3)
                               ##              mu1.est: estimated mean function in the 1st dimension (vector of length ngrid)  
                               ##              mu2.est: estimated mean function in the 2nd dimension (vector of length ngrid)  
                               ##              mu3.est: estimated mean function in the 3rd dimension (vector of length ngrid)  
                               ##              v.est: estimated eigen values (vector of length 3)
                               
){
  
#############################################################################
  ## Description: Function for the CS-MFPCA estimation algorithm Steps 3 in Section 4).
  ## Definition: ngrid: number of total time points;  
  ##             M: number of dimensions in the original multivariate functional outcome X(t);
  ##             N: number of total individuals.
  ## Args:        see above
  ## Returns:     list()
  ##              X0.est: estimates of transformed multivariate functional outcome in the 1st dimension (matrix of dimension N*ngrid)
  ##              X1.est: estimates of transformed multivariate functional outcome in the 2nd dimension (matrix of dimension N*ngrid)
  ##              X2.est: estimates of transformed multivariate functional outcome in the 3rd dimension (matrix of dimension N*ngrid)
  ##              X3.est: estimates of transformed multivariate functional outcome in the 4th dimension (matrix of dimension N*ngrid)
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
  
## Define trial points
gridPoints <- data.T$gridPoints


## Withdraw the estimated three-dimensional multivariate functional outcome
Z1.est <- MFPCAout$Z1.est
Z2.est <- MFPCAout$Z2.est
Z3.est <- MFPCAout$Z3.est

## Define number of individuals and number of total time points
N <- nrow(Z1.est)
ngrid <- ncol(Z1.est)

X0.est <- X1.est <- X2.est <- X3.est <- Z1.est
for (i in 1:N){
  for (j in 1:ngrid){
    X0.est[i,j] <- 1/(1+exp(Z1.est[i,j])+exp(Z2.est[i,j]) + exp(Z3.est[i,j]))
    X1.est[i,j] <- exp(Z1.est[i,j])/(1+exp(Z1.est[i,j])+exp(Z2.est[i,j]) + exp(Z3.est[i,j]))
    X2.est[i,j] <- exp(Z2.est[i,j])/(1+exp(Z1.est[i,j])+exp(Z2.est[i,j]) + exp(Z3.est[i,j]))
    X3.est[i,j] <- exp(Z3.est[i,j])/(1+exp(Z1.est[i,j])+exp(Z2.est[i,j]) + exp(Z3.est[i,j]))
  }
}


## Construct output
out <- list(X0.est = X0.est, X1.est = X1.est, X2.est = X2.est, X3.est = X3.est)
return(out)

}
  
