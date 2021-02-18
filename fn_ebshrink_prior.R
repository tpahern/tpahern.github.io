# Function to carry out empirical Bayes (EB) shrinkage

# Adapted from:
# Steenland K, Bray I, Greenland S, Boffetta P. Empirical Bayes adjustments for multiple results in hypothesis-generating or surveillance studies. Cancer Epidemiol Biomarkers Prev. 2000;9(9):895-903.

# bhat is a vector of estimates

# vhat is a vector of variances for the estimates in bhat

# niter is, obviously, the number of iterations to carry out (the default is 1,000)

# z is a matrix used to encode prior knowledge; if not using a prior, the default z is a column vector of length n

ebshrink_prior <- function(bhat, vhat, z=matrix(1,nrow=length(bhat),ncol=1), niter=1000){
  
  # define the number of estimates and number of params in prior model
  n <- length(bhat)
  p <- ncol(z)
  
  # set initial value of true variance to zero
  vart <- 0
  
  # create an empty vector to store estimates of true variance
  # set its first element to zero
  v <- rep(NA, niter+1)
  v[1] <- 0
  
  # generate a column of 1s for use in calculations
  I <- matrix(1, nrow=n, ncol=1)
  
  # start iterative loop
  for(i in 1:niter){
    
    # define weights and calculate their sum
    W <- solve(vhat + vart * diag(n))
    w <- t(I) %*% W %*% I
    
    # calculate overall mean (pi) and vector of prior means (mu)
    pi <- solve (t(z) %*% W %*% z) %*% t(z) %*% W %*% bhat
    mu <- z %*% pi
    
    # calculate D, Vhat_obs (varo) and Vhat_mean (varm)
    D <- bhat - mu
    varo <- t(D) %*% W %*% D/w
    varm <- (t(I) %*% W %*% vhat %*% I)/w
    
    # calculate new estimate of Var_true, set to zero if negative, and store in vector v
    vart <- (n/(n-p))*varo-varm
    vart <- max(vart,0)
    v[i+1] <- vart
    
    # stop iterating when Vhat_true has converged (at 10^-8) and calculate new estimates (newbeta) and variances (newvar)
    if((v[i+1] %% v) < 1e-007) {
      vart <- v[i+1]
      vstar <- (((n-p-2) * vhat)/(n-p))
      tstar <- vart * diag(n) + vhat - vstar
      newbeta <- W %*% ((tstar %*% bhat) + (vstar %*% mu))
    
      H <- z %*% (solve(t(z) %*% W %*% z)) %*% t(z) %*% W
      bstar <- ((n-p-2) * W %*% vhat)/(n-p)
      E <- bstar %*% D * sqrt(varm/diag(vhat))
      newvar <- diag(vhat-(t(diag(n)-H) %*% vhat %*% bstar) + (2*E %*% t(E))/(n-p-2))
      
      break
    }
  }
}