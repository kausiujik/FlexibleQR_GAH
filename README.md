# FlexibleBQR_GAH: Flexible Bayesian Quantile Regression based on the Generalized Asymmetric Huberised-type Distribution

This package implements Bayesian quantile regression based on the generalized asymmetric Huberised-type distribution as proposed by the following paper.

Hu, W. and Zhang W. (2024). 

Functions are implemented in *GAH_functions.R*.

## Simple Example

Install function and demo data.

```{r}
source("GAH_functions.R")
load("simdata.RData")
```

Five functions are included in *GAH_functions.R*.

- `g` defines the transformation function.
- `gammaINT` is used to calculate the finite end points for parameter $\gamma$. 
- `f_GAH` and `logGAH` are used to characterize the GAH distribution. 
- `BQR_GAH` fits the input data with the proposed method.

Input of `BQR_GAH`

- `y`: vector of continuous response variables
- `X`: (n, p)-matrix of covariates with intercept term (n: total sample size; p: the number of covariates)
- `tau0`: expected quantile level
- `nburn`: integer scalars indicate the numbers of burn-in samples
- `niter`: integer scalars indicate the numbers of posterior samples
- `thin`: integer scalars indicate the thinning rate
- `method`: selected from `alasso`(default) or `lasso`
- `beta0`: the initial value of $\pmb\beta$, `0` for default
- `gamma0`: the initial value of $\gamma$, `0` for default
- `eta0`: the initial value of $\eta$, `1` for default
- `rho20`: the initial value of $\rho^2$, `1` for default
- `a`: hyperparameter in Gamma prior for $\eta$, `1` for default
- `b`: hyperparameter in Gamma prior for $\eta$, `1` for default
- `r`: hyperparameter in Gamma prior for regularization, `1` for default
- `delta`: hyperparameter in Gamma prior for regularization, `1` for default
- `model.fit`: logical value indicate whether return DIC and LPML or not, `FALSE` for default

Output of `BQR_GAH`

- `beta`: posterior samples of $\pmb\beta$
- `gamma`: posterior samples of $\gamma$
- `eta`: posterior samples of $\eta$
- `rho2`: posterior samples of $\rho^2$
- `beta_brn`: burn-in samples of $\pmb\beta$
- `gamma_brn`: burn-in samples of $\gamma$
- `eta_brn`: burn-in samples of $\eta$
- `rho2_brn`: burn-in samples of $\rho^2$
- `beta_hat`(with `model.fit = TRUE`): posterior median of posterior samples of $\pmb\beta$
- `gamma_hat`(with `model.fit = TRUE`): posterior median of posterior samples of $\gamma$
- `eta_hat`(with `model.fit = TRUE`): posterior median of posterior samples of $\eta$
- `rho2_hat`(with `model.fit = TRUE`): posterior median of posterior samples of $\rho^2$
- `DIC`(with `model.fit = TRUE`): Deviance Information Criterion
- `LPML`(with `model.fit = TRUE`): Log-Pseudo Marginal Likelihood

```{r}
result_lasso <- BQR_GAH(y, X, tau0 = 0.25, nburn = 2000, niter = 10000, thin = 10, method = "lasso", model.fit = TRUE)
result_lasso$beta_hat    # regression coefficients
result_lasso$DIC         # DIC for the model
result_lasso$LPML        # LPML for the model

result_alasso <- BQR_GAH(y, X, tau0 = 0.25, nburn = 2000, niter = 10000, thin = 10, method = "alasso", model.fit = FALSE)
result_alasso$beta_brn   # burn-in samples of regression coefficients
result_alasso$beta       # posterior samples of regression coefficients
```
