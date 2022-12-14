# IRR
The goal of IRR is to assess interrater and intrarater reliability in the presence of nested data by fiting the Bayesian independent (BIN) model, Bayesian partially nested (BPN) model, and Bayesian fully nested (BFN) models.

We will be exploring the running gait (drone) dataset that was previously analyzed in Lafferty et al. (2022). The study consisted of 32 runners, evaluated by 3 
raters, over 2 timepoints, exploring 2 feet, in 2 locations yielding a total of 768 data points.  


## Installing packages and loading functions/data
```{r}
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(rstan)
library(dplyr)
library(knitr)
library(kableExtra)
library(loo)

source("compiling.R")
source("functions.R")
load("mydata.Rdata")
```

## Exploring the data
```{r}
knitr::kable(head(df.drone)) %>% kable_classic_2(full_width = F)
```
## Model Fitting

### Fitting BIN
```{r}
m1.BIN = model_BIN(df.drone, addFixedEff = TRUE, fixed_eff_intercept = TRUE, 
             beta_a = 5, beta_b = 5, gamma_a = 3, gamma_b = 1.5, beta_mean = 0, 
             beta_sigma = 1/0.3, niters = 1000, nwarmup = 200, nchains = 1)

print(m1.BIN,
      pars = c("beta_parm","sigma_S","sigma_R","sigma_T","Corr_R","Corr_T"))
loo(m1.BIN)
loo::waic(extract_log_lik(m1.BIN,parameter_name = "log_lik"))
```

### Fitting BPN
```{r}
m1.BPN = model_BPN(df.drone,cov_T_str = "unstructured", addFixedEff = TRUE, 
                   fixed_eff_intercept = TRUE, beta_a = 5, beta_b = 5, gamma_a = 3, 
                   gamma_b = 1.5, beta_mean = 0, beta_sigma = 1/0.3, betak_mean = 0, 
                   betak_sigma = 1/0.3, rho_T_eta = 1, rho_R_eta = 1, niters = 100, 
                   nwarmup = 20, nchains = 1) 
print(m1.BPN,
      pars = c("beta_parm","sigma_S","sigma_T","rho_T","Corr_T"))
loo(m1.BPN)
loo::waic(extract_log_lik(m1.BPN,parameter_name = "log_lik"))
```

### Fitting BFN
```{r}
m1.BFN = model_BFN(df.drone, cov_T_str = "unstructured", cov_R_str = "unstructured", 
                   addFixedEff = TRUE, fixed_eff_intercept = TRUE, beta_a = 5, 
                   beta_b = 5, gamma_a = 3, gamma_b = 1.5, beta_mean = 0, beta_sigma = 1/0.3, 
                   betak_mean = 0, betak_sigma = 1/0.3, rho_T_eta = 1, rho_R_eta = 1, 
                   niters = 100, nwarmup = 20, nchains = 1)  
print(m1.BFN,
      pars = c("beta_parm","sigma_S","sigma_T","rho_T","Corr_T"))
loo(m1.BFN)
loo::waic(extract_log_lik(m1.BFN,parameter_name = "log_lik"))
```


