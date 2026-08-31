# Capture Recapture - A Simulation Study

## :warning: NOTICE :warning:

This repository is under development and being converted to a joint CRAN library and analysis repository. This is not currently intended for public consumption.

[![DOI](https://zenodo.org/badge/822625163.svg)](https://doi.org/10.5281/zenodo.16580487)

## Overview

<img align="right" src="inst/extdata/crcsim.png" style="float" width="200">

To effectively tackle the overdose crisis, a nuanced understanding of Opioid Use Disorder (OUD) prevalence is crucial, both broadly and within targeted cohorts. Healthcare interactions provide estimates but may overlook those outside the healthcare system, leading to underestimation. Capture recapture (CRC) analysis is valuable in estimating prevalence by addressing underreporting in surveillance. Rich in history, the Syndemics Lab has utilized CRC procedures to estimate OUD prevalence in Massachusetts. Conventionally, a stepwise model selection process (MSP) is employed to identify the model that best fits the data. However, the MSP in estimating group-stratified prevalence is less explored, especially with sparse data. Additionally, the MSP works well with aggregate, uninformed data, however modern CRC methodologies are able to leverage relationships in the data to improve estimates and confidence interval bounds. This study uses simulations to investigate different MSPs for selecting conventional log-linear CRC models, with a focus on their ability to precisely estimate strata prevalence and the comparative impacts of estimators given covariates.

## Data

The data used by this project is all completely simulated data. As such, we have a helper function that is able to generate contingency tables without covariates suitable for the MSP case as well as row-level data containing categorical covariates for the more advanced estimators. While we utilize simulated data, because this package is intended for public consumption, one need only reshape their data into the expected format for their chosen estimator. Note that our `create_data` function does produce rows with "0" captures. These should be filtered out before running the model. We produce them for a true ground truth to compare against during our simulation runs, however these are the key estimated values via CRC.

Frequency Data:

| capture_1 | capture_2 | capture _3 | ... | N_ID |
|----------:|----------:|-----------:|----:|-----:|
| 0         | 0         | 1          | ... | 23   |
| 0         | 1         | 0          | ... | 42   |
| 0         | 1         | 1          | ... | 13   |
| 1         | 0         | 0          | ... | 56   |
| 1         | 0         | 1          | ... | 19   |
| 1         | 1         | 0          | ... | 74   |
| 1         | 1         | 1          | ... | 3    |

Covariate Data:

| capture_1 | capture_2 | capture _3 | --- | category_1 | category_2 | ... |
|----------:|----------:|-----------:|----:|-----------:|-----------:|----:|
| 0         | 0         | 1          | ... | 1          | 1          | ... |
| 1         | 1         | 0          | ... | 0          | 1          | ... |
| 0         | 1         | 1          | ... | 0          | 0          | ... |
| 1         | 0         | 0          | ... | 1          | 1          | ... |
| 0         | 0         | 1          | ... | 1          | 0          | ... |
| 1         | 1         | 0          | ... | 1          | 0          | ... |
| 1         | 0         | 1          | ... | 0          | 1          | ... |

## Creating Data

This library supports the ability to create data TODO: Fill in

## Generating Capture Probabilities

TODO: Fill in

## Generating Estimates

TODO: Fill in

## Analysis Notebooks

TODO: Fill in
