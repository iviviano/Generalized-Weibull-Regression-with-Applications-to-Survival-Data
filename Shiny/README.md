# A Shiny App for Weibull Regression on User Data

---

### Required Packages:
- [DT](https://cran.r-project.org/web/packages/DT/index.html)
- [mclm](https://cran.r-project.org/web/packages/mclm/readme/README.html)
- [rgenoud](https://cran.r-project.org/web/packages/rgenoud/index.html)
- [rjson](https://cran.r-project.org/web/packages/rjson/index.html)
- [rstan](https://mc-stan.org/users/interfaces/rstan)
- [shiny](https://www.rstudio.com/products/shiny/)
- [shinydashboard](https://rstudio.github.io/shinydashboard/index.html)
- [shinyvalidate](https://rstudio.github.io/shinyvalidate/)
- [survival](https://cran.r-project.org/web/packages/survival/index.html)
- [survminer](https://cran.r-project.org/web/packages/survminer/index.html)
- [tidyverse](https://www.tidyverse.org)

---

### Operation Details
The user is asked to input a data file. Several file types are supported. The user is asked to describe their data file. They must input:
- the number of continuous and discrete factors
- the column numbers of the data set corresponding to the failure times, censoring status, and each factor

There are several options for the model:
- Which model to run on the data
- The number of markov chains to be run
- The number of iterations to run each chain for
- Whether to use an orthogonal decomposition of the data matrix for computations
- Whether to use a normalization of the data for computations
None of these options should affect the parameter estimates besides possibly improving convergence. The user may down

Note: stan models sometimes have convergence issues. If you are getting very large values for the predictions from cross-validation (upwards of 10 orders of magnitude larger than the real data) try running cross-validation again, possibly with more iterations. 

---

### Implementation Details
- Bayesian models use NUTS and are implemented with the open source software [Stan](https://mc-stan.org/users/interfaces/rstan)
- Likelihood-based diagnostics use a population optimization algorithm implemented by [rgenoud](https://cran.r-project.org/web/packages/rgenoud/index.html)
- Discrete factors are handled by:
  1. Designating the first value as the baseline
  2. For each other observed value, an indicator covariate is formed of comparison against the baseline
- Cross Validation with censoring:
  1. Predictions are created for every data point
  2. The error statistics are calculated by omitting censored observations
  3. The prediction values used in error calculation are the mean of all the samples from the markov chain.
  4. The file of predictions contains both the matrix of samples and the mean values of the samples used for error calculations
- Cross Validation is not allowed if no factors are given

---

### Known Bugs
- Some of the survival plots are not showing up
- 
