# Generalized Weibull Regression with Applications to Survival Data
 ---
 
This is the github repository for the Bayesian Parameter Estimation REU project at the University of Michigan Dearborn in Summer of 2023. Abstract:

> We propose a Bayesian generalized Weibull regression method and develop Accelerated Time to Failure models using Bayesian methods. The parameter estimation procedure is carried out using Hamiltonian Monte Carlo algorithm with No-U-Turn Sampler and compares the results of generalized Weibull regression with exponentiated Weibull regression, Weibull regression, and log-normal distribution across simulated and clinical data sets. We examine the effectiveness of generalized Weibull distribution as a survival model and compare it to more studied probability distributions. In addition to monotone and bathtub hazard shapes, the additional shape parameter in the generalized Weibull distribution provides flexibility to model a broader class of monotone hazard rates. 
>
> Different model diagnostics tests such as the five-fold cross validation method, Akaike information criterion, and the Bayesian information criterion demonstrate the generalized Weibull distribution’s prediction accuracy. We found the generalized Weibull model to have superior results for multiple hazard function types. We also developed a Shiny app to dynamically visualize our models using real life data from lung cancer, heart attack, and German breast cancer studies. The app’s interactive interface lets users input their own data sets and see Weibull regression analysis on their data. 

<p markdown = "1">
The scripts used to simulate and analyze simulated data are located in the <a ref ="https://github.com/iviviano/Generalized-Weibull-Regression-with-Applications-to-Survival-Data/tree/main/Simulated-Data">Simulated Data</a> folder. All files required to run the shiny app are located in the <a href = "https://github.com/iviviano/Generalized-Weibull-Regression-with-Applications-to-Survival-Data/tree/main/Shiny">Shiny</a> folder.
</p>

---
# Acknowledgement

We would like to thank the National Science Foundation (DMS-1950102 and DMS-2243808), the National Security Agency (H98230-23), the College of Arts, Sciences, and Letters, and the University of Michigan Dearborn Department of Mathematics and Statistics for their support.
