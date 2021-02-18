# Function to carry out empirical Bayes (EB) shrinkage without a prior

# Adapted from:
# Steenland K, Bray I, Greenland S, Boffetta P. Empirical Bayes adjustments for multiple results in hypothesis-generating or surveillance studies. Cancer Epidemiol Biomarkers Prev. 2000;9(9):895-903.

# bhat is a vector of estimates

# vhat is a vector of variances for the estimates in bhat

# niter is, obviously, the number of iterations to carry out (the default is 100)

ebshrink_ig <- function(bhat, vhat, niter=100){

  # define the number of estimates
  n <- length(bhat)
  
  # initialize Vhat_true (vart) at zero
  vart <- 0
  
  # set initial value for Vhat_true (vart) to zero
  v <- rep(NA, niter + 1)
  v[1] <- 0
  
  # start iterative loop
  for(i in 1:niter){

    # define a vector of weights (w) and calculate mean parameter (pi)
    w <- 1/(vhat + vart)
    pi <- sum(w * bhat)/sum(w)
    
    # calculate D, Vhat_obs (varo), and Vhat_mean (varm)
    D <- bhat - pi
    varo <- sum(w * D * D)/sum(w)
    varm <- sum(w * vhat)/sum(w)
    
    # calculate new estimate of Var_true, set to 0 if neg, store in v
    vart <- max((varo - varm), 0)
    v[i+1] <- vart
    
    # terminate iterations when Vhat_true has converged
    # calculate new vectors of estimates and variances
    if((v[i+1] - v[i]) < 1e-7){
      variter <<- v
      vart <<- v[i+1]
      newbeta <<- w*((vart * bhat) + (vhat * pi))
      E <<- w * vhat * D * sqrt(varm/vhat)
      newvar <<- vhat * (1 - vhat * w) + (2 * E * t(E))/n
      break
    } #terminate
  } #1:niter
#  return(vart)
#  return(newbeta)
#  return(E)
#  return(newvar)
} #function

