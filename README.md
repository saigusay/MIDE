---
output:
  html_document: default
  pdf_document: default
---
# Robust Minimum Divergence Estimation in a Spatial Poisson Point Process

This repository contains all of the code needed to conduct the simulation studies and to reproduce the data analysis in the manuscript entitled "Robust minimum divergence estimation in a spatial Poisson point process".

## "code" Folder

This folder contains all of the R functions and scripts needed to implement the proposed estimation algorithm to generate the tables and figures included in the manuscript. 
The scripts require that the following package is installed in R: *MASS*. 
Once these packages are installed (the scripts will load them as needed), and once the working directory has been set to the location of the main repository on the local machine, the code should run without error. 

+ **functions.R**: Contains all functions needed to estimate the regression parameter for the spatial Poisson point process model. 
This script is sourced by all other simulation and data analysis scripts, and so should not need to be run separately.

+ **sample-simulation**: Conducts simulations for the settings of the true parameter values for target and contamination distribution.
This script generates the estimates of the regression parameters in a "result/XXX" subfolder, where the output folder is automatically created.
The result files are needed to produce the figures and tables in the manuscript.

+ **sample-simulation-cv**: Conducts simulations where the value of tuning parameter is selected based on the mean RTMSPE calculated through 2-fold cross-validation. 
This script generates the estimates of the regression parameters in a "result/XXX" subfolder, where the output folder is automatically created.
The result files are needed to produce the figures and tables in the manuscript.

