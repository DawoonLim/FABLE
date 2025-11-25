# FABLE (Factor Analysis with the BLEssing of dimensionality)

R package to implement Factor Analysis with BLEssing of dimensionality (FABLE) approach. 
Provides posterior samples of the covariance matrix modeled using factor analysis, along with providing posterior mean (without carrying out sampling).

All results in the paper can be replicated with scripts in `extras/replicationCodes` .

Paper: https://arxiv.org/abs/2404.03805# .

## Install 

Please install the `devtools` package in `R`. 

Use`devtools::install_github("shounakch/FABLE")` to install the package.

## Example Usage

```
set.seed(1)
  
Lambda = matrix(rnorm(p*k, mean = 0, sd = lambdasd), nrow = p, ncol = k)
BinMat = matrix(rbinom(p*k, 1, 1-pi0), nrow = p, ncol = k) 
Lambda = Lambda * BinMat
  
Sigma0 = runif(p, 0.5, 5)
  
M = matrix(rnorm(n*k), nrow = n, ncol = k)
E = matrix(rnorm(n*p), nrow = n, ncol = p)
E = sweep(E, 2, sqrt(Sigma0), "*")
  
Y = (M %*% t(Lambda)) + E
    
FABLEPostMean = FABLEPosteriorMean(Y, gamma0 = 1, delta0sq = 1, maxProp = 0.95)
FABLESamples = FABLEPosteriorSampler(Y, gamma0 = 1, delta0sq = 1, maxProp = 0.95, MC = 1000)
```
