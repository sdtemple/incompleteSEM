#' MLE for complete general stochastic epidemic model
#'
#' Compute MLE for complete epidemic observations.
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#'
#' @return MLE for (beta, gamma)
#'
#' @export
mle_complete_gsem = function(r, i, N){
  n = length(r)
  t = 0
  for(j in 1:n){
    t = t + sum(sapply(i, min, r[j]) - sapply(i, min, i[j]))
  }
  ri = sum(r - i)
  g = n / ri
  b = (n - 1) / (t + (N - n) * ri) * N
  return(c(b,g))
}

#' Likelihood for complete general stochastic epidemic
#' 
#' Compute exact likelihood for complete epidemic observations.
#' 
#' @param rates numeric vector: rates (beta, gamma)
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#' 
#' @return negative log likelihood
#' 
#' @export
likelihood_complete_gsem <- function(rates, r, i, N){
  
  beta <- rates[1]
  gamma <- rates[2]
  
  beta <- beta / N
  alpha <- which.min(i)
  n <- length(r)
  
  log_phi <- - beta * (N - n) * sum(r - i)
  log_f <- - gamma * sum(r - i) + n * log(gamma)
  
  psichi <- rep(0, n)
  for(j in (1:n)[-alpha]){
    rj <- r[j]
    ij <- i[j]
    X <- 0
    Y <- 0
    for(k in (1:n)[-j]){
      rk <- r[k]
      ik <- i[k]
      Y <- Y - beta * (min(rk, ij) - min(ik, ij))
      X <- X + as.numeric((ik < ij) & (ij < rk)) * beta
    }
    psichi[j] <- Y + log(X)
  }
  
  log_psichi <- sum(psichi[-alpha])
  #print(-log_phi)
  #print(-log_psichi)
  #print(-log_f)
  
  return(-(log_phi+log_f+log_psichi))
}
