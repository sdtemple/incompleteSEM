#' Simulate general stochastic epidemic model
#'
#' Draw infectious periods for general stochastic epidemic.
#'
#' @param rates numeric vector: rates (beta, gamma)
#' @param N integer: population size
#' @param m integer: positive shape
#'
#' @return matrix: infection times, removal times
#'
#' @export
rgsem = function(rates, N, m = 1){
  
  beta <- rates[1]
  gamma <- rates[2]
  
  # initialize vectors
  t = 0
  betaN = beta / N
  i = rep(Inf, N)
  r = rep(Inf, N)
  M = rep(0, N)
  alpha = sample(N, 1)
  i[alpha] = t
  
  # simulate epidemic
  St = sum(is.infinite(i))
  It = sum(is.finite(i)) - sum(is.finite(r))
  while(It > 0){
    
    # simulate time
    irate = betaN * It * St
    rrate = gamma * It
    t = t + rexp(1, rate = irate + rrate)
    
    # simulate transition
    x = rbinom(1, size = 1, prob = irate / (irate + rrate))
    if(x){
      # infect a susceptible
      if(St > 1){
        argx = sample(which(is.infinite(i), arr.ind = T), 1)
      } else{
        argx = which(is.infinite(i), arr.ind = T)
      }
      i[argx] = t
    } else{
      # remove an infected
      if(It > 1){
        argx = sample(which(is.infinite(r) & is.finite(i), arr.ind = T), 1)
      } else{
        argx = which(is.infinite(r) & is.finite(i))
      }
      M[argx] = M[argx] + 1
      if(M[argx] == m){ # after m renewals
        r[argx] = t
      }
    }
    
    # update (S,I) counts
    St = sum(is.infinite(i))
    It = sum(is.finite(i)) - sum(is.finite(r))
    
  }
  
  # formatting
  output = matrix(c(i,r),
                  nrow = N,
                  ncol = 2,
                  byrow = F)
  colnames(output) = c('i','r')
  return(output)
}

#' Sort epidemic data
#' 
#' Sort matrix by increasing removal times.
#' 
#' @param epi matrix: infection times, removal times
#' 
#' @return matrix: infection times, removal times
#' 
#' @export
sort_gsem <- function(epi){
  r <- epi[,2]
  i <- epi[,1]
  ind <- order(r)
  r <- r[ind]
  i <- i[ind]
  return(cbind(i,r))
}

#' Filter epidemic data
#' 
#' Keep cases only.
#' 
#' @param epi matrix: infection times, removal times
#' 
#' @return matrix: infection times, removal times
#' 
#' @export
filter_gsem <- function(epi){
  r <- epi[,2][is.finite(epi[,2])]
  i <- epi[,1][is.finite(epi[,2])]
  return(cbind(i,r))
}

#' Missing epidemic data
#' 
#' Impute NAs for infection and removal times
#' 
#' @param epi infection and removal times
#' @param p expected proportion of complete pairs observed
#' @param q probability infection time missing
#' 
#' @return matrix: infection times, removal times
#' 
#' @export 
decomplete_gsem <- function(epi, p, q = 1){
  r <- epi[,2]
  i <- epi[,1]
  n <- length(r)
  # alpha <- which.min(i)
  # i[alpha] <- NA
  # for(j in (1:n)[-alpha]){
  for(j in 1:n){
    if(rbinom(1, 1, 1 - p)){
      if(rbinom(1, 1, q)){
        i[j] <- NA
      } else{
        r[j] <- NA
      }
    }
  }
  return(cbind(i,r))
}