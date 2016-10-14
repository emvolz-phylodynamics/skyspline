# skyspline

Phylodynamic inference using birth-rate splines. More documentation and vignettes coming soon.

This package includes functions for MLE/bootstrap and Bayesian MCMC estimation of population size and birth rate through time given genealogical data. 

Use fit.skyspline.ml for maximum likelihood or maximum a posterior inference. 
parboot.skyspline.mle can be used to estimate CIs. 
skyspline.metrop.hastings provides a simple Metropolis-Hastings sampler for Bayesian inference. 
