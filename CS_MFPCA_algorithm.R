CS_MFPCA_algorithm <- function(data.obs # Output from function CS_MFPCA_simulation 
                                    # list of observed values, including:
                                     #           X0: original multivariate functional outcome in the reference dimension (matrix of dimension N*ngrid)
                                     #           X1: original multivariate functional outcome in the 1st dimension (matrix of dimension N*ngrid)
                                     #           X2: original multivariate functional outcome in the 2nd dimension (matrix of dimension N*ngrid) 
                                     #           X3: original multivariate functional outcome in the 3rd dimension (matrix of dimension N*ngrid)
                                     #           Z1: transformed multivariate functional outcome in the 1st dimension (matrix of dimension N*ngrid)
                                     #           Z2: transformed multivariate functional outcome in the 2nd dimension (matrix of dimension N*ngrid) 
                                     #           Z3: transformed multivariate functional outcome in the 3rd dimension (matrix of dimension N*ngrid) 
                               
){
  
#############################################################################
## Description: Function for the CS-MFPCA estimation algorithm steps 2 in Section 4).
## Definition: ngrid: number of total time points;  
##             M: number of dimensions in the original multivariate functional outcome X(t);
##             N: number of total individuals.
## Args:        see above
## Returns:     list()
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

## Withdraw generated four-dimensional multivariate functional outcome X(t)
X0 <- data.obs$X0
X1 <- data.obs$X1
X2 <- data.obs$X2
X3 <- data.obs$X3


#################################
# CS-MFPCA ######################
#################################
## Create the transformed three-dimensional multivariate functional outcome Z(t)
Z1 <- log(X1 / X0)
Z2 <- log(X2 / X0)
Z3 <- log(X3 / X0)

# Create `funData` objects in each dimension of three-dimensional multivariate functional outcome Z(t)
Z1.obj <- funData(argvals = gridPoints, X = Z1) 
Z2.obj <- funData(argvals = gridPoints, X = Z2) 
Z3.obj <- funData(argvals = gridPoints, X = Z3) 

# Create `MultiFunData` objects from all dimensions of multivariate functional outcome Z(t)
Z.mv.obj <- multiFunData(Z1.obj, Z2.obj, Z3.obj)

# MFPCA upon the `MultiFunData` objects from all dimensions of multivariate functional outcome Z(t)
uniExpansions <- list(list(type = "uFPCA"), 
                      list(type = "uFPCA"),
                      list(type = "uFPCA"))
CS.MFPCA.out <- MFPCA(Z.mv.obj, M = 3, uniExpansions = uniExpansions, fit = TRUE)

# Withdraw the estimated transformed multivariate functional outcome Z(t) from MFPCA
Z1.est <- CS.MFPCA.out$fit[[1]]@X
Z2.est <- CS.MFPCA.out$fit[[2]]@X
Z3.est <- CS.MFPCA.out$fit[[3]]@X

# Withdraw the estimated multivariate eigenfunctions for Z(t) from MFPCA
phi1.est <- t(CS.MFPCA.out$functions[[1]]@X)
phi2.est <- t(CS.MFPCA.out$functions[[2]]@X)
phi3.est <- t(CS.MFPCA.out$functions[[3]]@X)


# Withdraw the estimated mean functions for Z(t) from MFPCA
mu1.est <- CS.MFPCA.out$meanFunction[[1]]@X
mu2.est <- CS.MFPCA.out$meanFunction[[2]]@X
mu3.est <- CS.MFPCA.out$meanFunction[[3]]@X

# Withdraw the estimated eigen values for Z(t) from MFPCA
v.est <- CS.MFPCA.out$values

## Construct output
out <- list(Z1.est = Z1.est, Z2.est = Z2.est, Z3.est = Z3.est,
            phi1.est = phi1.est, phi2.est = phi2.est, phi3.est = phi3.est,
            mu1.est = mu1.est, mu2.est = mu2.est, mu3.est = mu3.est,
            v.est = v.est)

return(out)
}



