#' PBLA for partial general stochastic epidemic model
#'
#' Compute approximate likelihood for partial epidemic observations.
#'
#' @param rates numeric vector: rates (beta, gamma)
#' @param r numeric vector: removal times
#' @param N integer: population size
#' @param m integer: positive shape
#' @param A integer: plausible patient zeros
#' @param lag numeric: fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_partial_gsem = function(rates, r,  N, m = 1, A = 1, lag = 0){
  
  beta <- rates[1]
  gamma <- rates[2]
  
  # initialize
  n = length(r)
  r1 = r[1]
  beta = beta / N
  
  # change of variable to delta
  if(n < N){
    B = beta * (N - n)
    delta = gamma + B
  } else{ # handles entire population infected
    if(n == N){delta = gamma}
  }
  
  # calculate log likelihood (line six)
  ia = rep(-log(A), A)
  ip = - delta * (r[1:A] - r1)
  z = ia + ip
  
  # copy and paste from is.integer documentation
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }
  
  if((any(beta < 0)) | (any(gamma < 0)) |
     (!is.wholenumber(N)) | (N <= 0) |
     (!is.wholenumber(m)) | (m <= 0) |
     (!is.wholenumber(A)) | (A <= 0)){
    # invalid parameters
    return(1e15)
  } else{
    if(m == 1){ # exponential infectious periods
      
      # evaluate psi and chi terms
      psichi = rep(0, n)
      b = beta
      denom = 2 * delta * (b + delta)
      for(j in (1:n)){
        X = 0
        Y = 0
        rj= r[j]
        for(k in (1:n)[-j]){
          rk = r[k]
          # lemma 1
          if(rj - lag < rk){
            w = delta / denom  * exp(- delta * (rk - (rj - lag)))
            x = delta * w
            y = 1 - b * w
          } else{
            w = delta / denom * exp(- delta * ((rj - lag) - rk))
            x = delta * w
            y = delta / (b + delta) + b * w
          }
          # line twelve
          X = X + b * x / y
          Y = Y + log(y)
        }
        #print(c(j, log(X), Y))
        psichi[j] = Y + log(X)
      }
      
      # line eight
      for(alpha in 1:A){z[alpha] = z[alpha] + sum(psichi[-alpha])}
      z = matrixStats::logSumExp(z)
      a = n * log(gamma / delta)
      
      # negative log likelihoods
      return(-(a+z))
    } else{ # erlang case
      
      # evaluate psi and chi terms
      psichi = rep(0, n)
      b = beta
      for(j in (1:n)){
        X = 0
        Y = 0
        rj= r[j]
        for(k in (1:n)[-j]){
          rk = r[k]
          # lemma 4
          if(rj - lag < rk){
            U = 0
            V = 0
            for(l in 0:(m-1)){
              v = 0
              for(p in 0:l){
                v = v + choose(m + p - 1, p) /
                  factorial(l - p) *
                  ((rk - rj + lag) ^ (l - p)) /
                  ((delta + delta) ^ (m + p))
              }
              U = U + v / ((delta + b) ^ (m - l))
              V = V + v * (delta ^ l) *
                (((delta / (delta + b)) ^ (m - l)) - 1)
            }
            w = exp(- delta * (rk - rj + lag)) * (delta ^ m)
            x = (delta ^ m) * w * U
            y = 1 + w * V
          } else{
            U = 0
            V = 0
            for(l in 0:(m-1)){
              v = 0
              for(p in 0:(m-1)){
                v = v + choose(l + p, p) /
                  factorial(m - p - 1) *
                  ((rj - lag - rk) ^ (m - p - 1)) /
                  ((delta + delta) ^ (l + p + 1))
              }
              U = U + v / ((delta + b) ^ (m - l))
              V = V + v * (delta ^ l) *
                (((delta / (delta + b)) ^ (m - l)) - 1)
            }
            w = exp(- delta * (rj - lag - rk)) * (delta ^ m)
            x = (delta ^ m) * w * U
            y = 1 + (w * V) -
              pgamma(rj - lag - rk, m, delta) *
              (1 - ((delta / (delta + b)) ^ m))
          }
          # line twelve
          X = X + b * x / y
          Y = Y + log(y)
        }
        psichi[j] = Y + log(X)
      }
      
      # line eight
      for(alpha in 1:A){z[alpha] = z[alpha] + sum(psichi[-alpha])}
      z = matrixStats::logSumExp(z)
      a = n * m * log(gamma / delta)
      # negative log likelihoods
      return(-(a+z))
    }
  }
}

#' PBLA for stochastic epidemic model
#'
#' Compute pair-based likelihood approximation. Supports Erlang infectious periods.
#'
#' @param r numeric vector: removal times
#' @param beta matrix of rates
#' @param gamma numeric vector of rates
#' @param m positive integer shape
#' @param A integer patient zeros
#' @param lag numeric fixed lag
#'
#' @return negative log likelihood
#'
#' @export
pbla_partial = function(r, beta, gamma, m = 1, A = 1, lag = 0){
  
  # initialize
  n = length(r)
  N = ncol(beta)
  r1 = r[1]
  
  # change of variable to delta
  if((n < (N - 1)) & (n > 1)){
    B = apply(beta[1:n,(n+1):N], 1, sum)
    delta = gamma + B
  } else{ # handles special cases
    if(n == N){delta = gamma}
    if(n == (N - 1)){delta = gamma + beta[1:(N-1),N]}
    if(n == 1){delta = gamma + sum(beta[1,2:N])}
  }
  
  # calculate log likelihood (line six)
  ia = rep(-log(A), A)
  ip = - delta[1:A] * (r[1:A] - r1)
  z = ia + ip
  
  # copy and paste from is.integer documentation
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }
  
  if((any(beta < 0)) | (any(gamma < 0)) |
     (!is.wholenumber(m)) | (m <= 0) |
     (!is.wholenumber(A)) | (A <= 0)){
    # invalid parameters
    return(1e15)
  } else{
    if(m == 1){ # exponential infectious periods
      
      # evaluate psi and chi terms
      psichi = rep(0, n)
      for(j in (1:n)){
        X = 0
        Y = 0
        rj= r[j]
        deltaj = delta[j]
        for(k in (1:n)[-j]){
          b = beta[k,j]
          rk = r[k]
          deltak = delta[k]
          denom = (deltaj + deltak) * (b + deltak)
          # lemma 1
          if(rj - lag < rk){
            w = deltaj / denom  * exp(- deltak * (rk - (rj - lag)))
            x = deltak * w
            y = 1 - b * w
          } else{
            w = deltak / denom * exp(- deltaj * ((rj - lag) - rk))
            x = deltaj * w
            y = deltak / (b + deltak) + b * w
          }
          # line twelve
          X = X + b * x / y
          Y = Y + log(y)
        }
        psichi[j] = Y + log(X)
      }
      
      # line eight
      for(alpha in 1:A){z[alpha] = z[alpha] + sum(psichi[-alpha])}
      z = matrixStats::logSumExp(z)
      a = sum(log(gamma / delta))
      
      # negative log likelihoods
      return(-(a+z))
    } else{ # erlang case
      
      # evaluate psi and chi terms
      psichi = rep(0, n)
      for(j in (1:n)){
        X = 0
        Y = 0
        rj= r[j]
        deltaj = delta[j]
        for(k in (1:n)[-j]){
          b = beta[k,j]
          rk = r[k]
          deltak = delta[k]
          # lemma 4
          if(rj - lag < rk){
            U = 0
            V = 0
            for(l in 0:(m-1)){
              v = 0
              for(p in 0:l){
                v = v + choose(m + p - 1, p) /
                  factorial(l - p) *
                  ((rk - rj + lag) ^ (l - p)) /
                  ((deltaj + deltak) ^ (m + p))
              }
              U = U + v / ((deltak + b) ^ (m - l))
              V = V + v * (deltak ^ l) *
                (((deltak / (deltak + b)) ^ (m - l)) - 1)
            }
            w = exp(- deltak * (rk - rj + lag)) * (deltaj ^ m)
            x = (deltak ^ m) * w * U
            y = 1 + w * V
          } else{
            U = 0
            V = 0
            for(l in 0:(m-1)){
              v = 0
              for(p in 0:(m-1)){
                v = v + choose(l + p, p) /
                  factorial(m - p - 1) *
                  ((rj - lag - rk) ^ (m - p - 1)) /
                  ((deltaj + deltak) ^ (l + p + 1))
              }
              U = U + v / ((deltak + b) ^ (m - l))
              V = V + v * (deltak ^ l) *
                (((deltak / (deltak + b)) ^ (m - l)) - 1)
            }
            w = exp(- deltaj * (rj - lag - rk)) * (deltaj ^ m)
            x = (deltak ^ m) * w * U
            y = 1 + (w * V) -
              pgamma(rj - lag - rk, m, deltaj) *
              (1 - ((deltak / (deltak + b)) ^ m))
          }
          # line twelve
          X = X + b * x / y
          Y = Y + log(y)
        }
        psichi[j] = Y + log(X)
      }
      
      # line eight
      for(alpha in 1:A){z[alpha] = z[alpha] + sum(psichi[-alpha])}
      z = matrixStats::logSumExp(z)
      a = sum(m * log(gamma / delta))
      # negative log likelihoods
      return(-(a+z))
    }
  }
}
