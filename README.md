# jsPGA
Predicting fish responses to climate change using a joint species, spatially dependent physiologically guided abundance model (https://doi.org/10.1002/ecy.4362)

## Information
This repository is cloned from the official USGS repository.

Repository Type: R and stan scripts 

Year of Origin:  2023

Year of Version: 2023

Version: 1.0.0

Digital Object Identifier (DOI): https://doi.org/10.5066/P959EMT5

USGS Information Product Data System (IPDS) no.: IP-158441 (internal agency tracking)

## Suggested Citation for Software

Custer, C., North, J. S., Schliep, E. M., and Wagner, T. 2023. Predicting fish responses to climate change using a joint species, spatially dependent physiologically guided abundance model. U.S. Geological Survey software release, https://doi.org/10.5066/P959EMT5


## Abstract

Predicting the effects of warming temperatures on the abundance and distribution of organisms under future climate scenarios often requires extrapolating species-environment correlations to thermal conditions not currently experienced by a species or to temperatures that exceed the range of observed data. For poikilotherms, incorporating species’ thermal physiology to inform extrapolations under novel thermal conditions can result in more realistic predictions.  Here, we  present a joint species, spatially dependent physiologically guided abundance (jsPGA) model for predicting multispecies responses to climate warming. The jsPGA uses a basis function approach to capture both species and spatial dependencies.  We first show that the jsPGA model accurately estimates parameters through a simulation study.  We then apply the jsPGA to predict the response of eight fish species observed across thousands of lakes in Minnesota, USA to projected climate warming. The jsPGA provides a new tool for predicting changes in abundance, distribution, and extinction probability of poikilotherms under novel thermal conditions.

This repository provides all the necessary R and Stan scripts to fit the jsPGA model

File descriptions:

R/Data Manipulation.R provides code to prepare the raw data for modeling.  Creates a few versions of a MNfish.rds object. Comments within the script provide details on differences between produced data objects.

R/sim_study.R reads in the standardized `MNfishz.rds` object and produces a simulated dataset that utilizes the observed data and predictor variables (i.e., environmental and effort) with simulated parameters to created simulated catch values.  
This file creates a list `simdat_M--_TC.rds` (where -- is the number of basis vectors specified within the code) that is suitable to use within the `cmdstanr` R package.

R/MN_case_study.R reads in the standardized `MNfishz.rds` object and produces a list `MNdat_M--_TC.rds` (where -- is the number of basis vectors specified within the code) that is suitable to use within the `cmdstanr` R package.

R/sim_analysis.R first merges all posteriors obtained from each model fit realization into a single posterior for inference.  
It then provides the code used to produce all the posterior summaries, tables, and figures for the simulated models found within the jsPGA paper.

R/MN_summary.R first merges all posteriors obtained from each model fit realization into a single posterior for analysis.  
It then provides the code used to produce the summaries, tables, and figures for the parameter estimates from the Minnesota case study found within the jsPGA paper.

R/MN_predictions.R produces the predictions for the thermal performance curves and relative abundance and their respective summaries, tables, and figures found within the jsPGA paper.

R/jsPGA_cmdstan.R is the code the reads in the list `.rds` objects from R/sim_study.R and R/MN_case_study.R and fits the jsPGA model within the Stan framework using the `cmdstanr` R package.

Stan/jspga.stan is the Stan code for the jsPGA model.
