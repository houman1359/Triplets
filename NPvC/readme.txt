glmnet: package to compute the GLM fits.

Ker_Copula: computes the density function using vine copulas for a given multivariate vector. An example of how the data should be prepared and fits is given in the file Simulation_Pairwise_git.m in the folder SIMULATIONS.
The main files are

Fit_vCopula.m      - Fit and evaluate kernel vine-Copula to the data

prep_vine.m        - prepares data to structure according to the nodes of vine and be used for the fits and density estimations.

predict_response.m - Predict mle and em estimates for the response.


