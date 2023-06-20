CONTENTS OF THIS FOLDER ——————————————

CS_MFPCA_tutorial.R : A step-by-step implementation of CS-MFPCA and the associated procedures described in "Constrained multivariate functional principal components analysis for novel outcomes in eye-tracking experiments".

CS_MFPCA_simulation.R : Function for simulating one data set under the simulation design described in Section 6 of "Constrained multivariate functional principal components analysis for novel outcomes in eye-tracking experiments".

CS_MFPCA_algorithm.R : Function for fitting multivariate functional principal components analysis (MFPCA) on the transformed multivariate functional outcome Z(t) included in Step 2 of the CS-MFPCA algorithm described in "Constrained multivariate functional principal components analysis for novel outcomes in eye-tracking experiments".

CS_MFPCA_prediction.R : Function for obtaining the predicted original multivariate functional outcome X(t) included in Step 3 of the CS-MFPCA algorithm described in "Constrained multivariate functional principal components analysis for novel outcomes in eye-tracking experiments".

INTRODUCTION ——————————————

The contents of this folder allow for implementation of the CS-MFPCA described in "Constrained multivariate functional principal components analysis for novel outcomes in eye-tracking experiments". Users can simulate a sample data set of the original four-dimensional multivariate functional outcome (CS_MFPCA_simulation.R) and apply the proposed CS-MFPCA algorithm (CS_MFPCA_algorithm.R). Also, we include tools to perform predictions for multvariate trajectories (CS_MFPCA_prediction.R). Detailed instructions on how to perform the aforementioned CS-MFPCA procedures, predictions of multivariate trajectories, and visualize results are included in CS_MFPCA_tutorial.R.

REQUIREMENTS ——————————————

The included R programs require R 4.0.5 (R Core Team, 2021) and the packages listed in CS_MFPCA_tutorial.R.

INSTALLATION ——————————————

Load the R program files into the global environment and install required packages using commands in CS_MFPCA_tutorial.R.
