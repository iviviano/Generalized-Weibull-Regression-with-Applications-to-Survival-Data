# Simulated Data

---

### File Contents

- SimulateData.R is an R script to simulate 4 data sets from the generalized gamma distribution. It applies a right censoring mechanism with various censoring rates to each data sets, and saves the resulting data in sim_data.csv
- sim_data.csv is a .csv file containing the simulated data
- AnalyzeSimulatedData.Rmd contains the cross-validation functions used to compare the prediction accuracies of different models

---

### Implementation Details

The generalized gamma distribution is sampled using the [flexsurv](https://cran.r-project.org/web/packages/flexsurv/index.html) package in R. Bathtub hazard data is simulated with parameter values: ...

---

### Handling Censoring

The following algorithm is used to apply right censoring to the simulated data sets:

1. Take a random sample of the data points to censor
2. For each data point, choose a uniform random value between 0 and the true failure time. This becomes the censored value.
3. Mark censored observations as such

A function implementing this algorithm is contained in SimulateData.R.