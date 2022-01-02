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
