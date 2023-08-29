
# **Appendix B: Hierarchical model code and details on Bayesian Lasso**

 

## Hierarchical model code

``` r
# This model estimates point-level occupancy probability as a joint 
# probability of block-level human footprint and point-level habitat, 
# where blocks contain >= 1 point. 

sink("block_occupancy_lasso.txt")
cat("
  model{
  #Priors on block-level coefficients (for Bayesian Lasso, uses double exponential distribution)
  beta.block[1] ~ dnorm(0, 1)
  for(i in 2:n.beta.block){
    beta.block[i] ~ ddexp(0, lambda) # where lamdba is the double exponential variance parameter specified for each Bayesian Lasso iteration
  }

  #Priors on point-level coefficients
  beta.point[1] ~ dnorm(0, 1)
  for(i in 2:n.beta.point){
    beta.point[i] ~ dnorm(0, 0.1)
  }
  
  # Priors on detection coefficients 
  for(i in 1:n.beta.p){
    beta.p[i] ~ dnorm(0, 0.1)
  }
  
  # Likelihood on mean within-block occupancy probability
  for(i in 1:n.block){
    logit(psi[i]) <- inprod(beta.block, cov.block[i, ]) # block-level occupancy probability
    delta[i] <- 1 - pow((1 - psi[i]), (1/n.scale)) # the proportion of area occupied as a function of human footprint / constant point-level occupancy probability
    N[i] ~ dbinom(delta[i], n.scale) # Number of occupied points as a function of block-level human footprint
  }
  
  # Likelihood point-level occupancy as a function of habitat, conditional on block-level occupancy rate
  for(i in 1:n.site){
    zb[i] ~ dbern(delta[point.block[i]]) # Point-level occupancy as a function of block-level parameters
    logit(theta[i]) <- inprod(beta.point, cov.point[i, ])
    zp[i] ~ dbern(zb[i]*theta[i])
    # Likelihood on detecting the species conditional on occupancy
    for(j in 1:n.surv[i]){
      logit(p[i, j]) <- inprod(beta.p, cov.p[i, j, ])
      y[i, j] ~ dbern(zp[i]*p[i, j])
    
      prob1[i, j] <- pow(p[i, j], y[i, j])*pow((1-p[i, j]), (1 - y[i, j]))
      prob2[i, j] <- 1 - p[i, j]
    }
    # Model likelihood calculation 
    term1[i] <- indicator[i]*(psi[point.block[i]]*theta[i])*prod(prob1[i, 1:n.surv[i]])
    term2[i] <- (1 - indicator[i])*((1 - (psi[point.block[i]]*theta[i])) + (psi[point.block[i]]*theta[i])*prod(prob2[i, 1:n.surv[i]]))
    prob.y[i] <- term1[i] + term2[i]
    lprob.y[i] <- log(prob.y[i])
    
    # Calculate the probability of detecting it at least once, 
    prob.i[i] <- (zp[i]*theta[i])*(1 - prod(1 - p[i, 1:n.surv[i]]))  # Probability of occupancy * probability of detecting at least once
    zc[i] ~ dbern(prob.i[i])
  }
  
  #Derived quantities
  l.score <- -2*sum(lprob.y[]) # Model likelihood estimate
  e.count.total <- sum(zc[]) # Total number of occupied points 
}
", fill = TRUE)
```

 

## **Details on Bayesian Lasso implementation**

# We implemented the Bayesian Lasso procedure for each model type (additive, interactive, and total human footprint) at each scale using all coefficients and interactions on the block-level portion of the model, with specific variables depending on which model type was being estimated. Similar to Gerber et al. \[-@gerber2015\] and Stevens and Conway \[-@Stevens2019\], we searching over 50 potential values of prior variance for the Laplace distribution ranging from 0.1–5 on the log scale, which translated to ~1.1–148 when transformed. We selected the optimal prior for each model type at each scale by comparing models using Watanabe-Akaike Information Criterion \[WAIC, @watanabe2010\].

 

# For implementation in JAGS, we specified the double exponential distribution as the prior distribution for all block-level variables as:

# ‘beta.block\[i\] ~ ddexp(0, lambda)’ (see model code above), where lambda was the prior variance value for each realization of the model.

 

# We generated the lamdba values as:

log.lambda=seq(0.1, 5, length = 50)  
lambda \<- exp(log.lambda)  
lam \<- as.numeric() \# Which lamdba value to choose (1-50)  
mod.lam \<- lambda\[lam\]

 