# Title: Variable-wise kernel fuzzy c-means clustering algorithms with kernalization of the metric
# Author: Gabriel Teotonio
# Date: 2021/01/29

# Aux. functions -----
kernel_gaussian <- function(x, v) {
  
  dif <- (x - v)
  percentiles <- quantile(abs(dif)^2, names = FALSE, probs = c(.1, .9))
  two_sigma <- sum(percentiles)/2
  
  res <- c(length(x))
  for (k in 1:length(x)) {
    res[k] <- exp(-dif[k]^2)/two_sigma
  }
  
  return(res)
}

metric_membership <- function(u, m, x, v) {
  
  kernel_value <- kernel_gaussian(x, v)
  res <- sum(u^m * (2 * (1 - kernel_value)))
  
  return(res)
}



# Algorithm -----
vkfcm_k_lp <- function(x, c = 2, m = 1.1, T_limit = 150, error_cond = 10^(-10)) {
  
  # Initiate objects -----
  V <- matrix(0L, nrow = c, ncol = ncol(x))  
  U <- diag(1, nrow = c, ncol = nrow(x))
  L <- matrix(1, nrow = c, ncol = ncol(x))
  aux_prod <- c(1, 1)
  
  initials_kick <- sample(1:nrow(x), c, replace = FALSE)
  V <- x[initials_kick, ]
  
  while(error > error_cond || t <= T_limit) {
      
    # Update prototype -----
    for (i in 1:c) {
     for(j in 1:ncol(x)) {
       V[i, j] <- sum(U[i, ]^m * kernel_gaussian(x[, j], v[i, j]) * x[, j])/(U[i, ]^m * kernel_gaussian(x[, j], v[i, j]))
       
       aux_prod[i] <- aux_prod[i] * metric_membership(U[i, ], m, x[, j], v[i, j])
     }
    }
    
    # Update  weight ----
    for(i in 1:c) {
     for(j in 1:ncol(x)) {
       L[i, j] <- (aux_prod[i])^(1/ncol(x))/metric_membership(U[i, ], m, x[, j], v[i, j])
     }
    }
    
    # Update membership ----
    for(i in 1:c) {
     for(k in 1:nrow(x)) {
       ()^(-1)
     }
    }
  }
}



