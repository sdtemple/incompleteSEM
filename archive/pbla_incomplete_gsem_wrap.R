#' Optimization Wrapper for `pbla_incomplete_gsem`
#'
#' Run a pair-based likelihood approximation for a general stochastic epidemic model. Compatible with `pbla_uni`.
#'
#' @param rates numeric vector of rates
#' @param r numeric vector of (increasing) removal times
#' @param i numeric vector of (paired) infection times
#' @param N integer population size
#' @param etc other parameters to pass (e.g. A)
#'
#' @return negative log likelihood
#'
#' @export
pbla_uni_wrap <- function(rates, r, i, N, etc = NULL){
  beta = rates[1]
  gamma = rates[2]
  if(is.null(etc)){ # use defaults
    return(do.call(pbla_incomplete, list(r=r,i=i,beta=beta,gamma=gamma,N=N)))
  } else{ # pass in all parameters
    return(do.call(pbla_incomplete, c(list(r=r,i=i,beta=beta,gamma=gamma,N=N), etc)))
  }
}

pbla_two_step <- function(beta, r, i, N, etc = NULL){
  gamma <- mle_removal_rate(r, i)[1]
  if(is.null(etc)){ # use defaults
    return(do.call(pbla_uni_wrap, list(rates=c(beta,gamma),r=r,i=i,N=N)))
  } else{ # pass in all parameters
    return(do.call(pbla_uni_wrap, c(list(rates=c(beta,gamma),r=r,i=i,N=N), etc)))
  }
}
