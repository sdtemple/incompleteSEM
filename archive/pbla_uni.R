#' Universal PBLA
#' 
#' Compute pair-based likelihood approximation for 9 cases. Assumes exponential infectious periods.
#' 
#' @param r numeric vector of (increasing) removal times
#' @param i numeric vector of (paired) infection times
#' @param beta numeric rate
#' @param gamma numeric rate
#' @param N integer population size
#' @param A integer (possible) patient zeros 
pbla_uni <- function(r, i, beta, gamma, N, A = 1){
  
  ##### internal functions #####
  
  E_psi <- function(rk, rj, ik, ij, deltak, deltaj, beta){
    if(is.na(rk)){
      if(is.na(rj)){ # ik, ij
        return(E_psi_ik_ij(ik, ij, deltak, deltaj, beta))
      } else if(is.na(ij)){ # rj, ik
        return(E_psi_rj_ik(rj, ik, deltak, deltaj, beta))
      } else{ # rj, ik, ij
        return(E_psi_rj_ik_ij(rj, ik, ij, deltak, deltaj, beta))
      }
    }
    
    if(is.na(rj)){
      if(is.na(ik)){ # rk, ij
        return(E_psi_rk_ij(rk, ij, deltak, deltaj, beta))
      } else{ # rk, ik, ij
        return(E_psi_rk_ik_ij(rk, ik, ij, deltak, deltaj, beta))
      }
    }
    
    if(is.na(ij)){
      if(is.na(ik)){ # rk, rj
        return(E_psi_rk_rj(rk, rj, deltak, deltaj, beta))
      } else{ # rk, rj, ik
        return(E_psi_rk_rj_ik(rk, rj, ik, deltak, deltaj, beta))
      }
    }
    
    if(is.na(ik)){ # rk, rj, ij
      return(E_psi_rk_rj_ij(rk, rj, ij, deltak, deltaj, beta))
    }
    
    return(E_psi_rk_rj_ik_ij(rk, rj, ik, ij, deltak, deltaj, beta)) # rk, rj, ik, ij
  }
  
  E_psi_ind <- function(rk, rj, ik, ij, deltak, deltaj, beta){
    if(is.na(rk)){
      if(is.na(rj)){ # ik, ij
        return(E_psi_ind_ik_ij(ik, ij, deltak, deltaj, beta))
      } else if(is.na(ij)){ # rj, ik
        return(E_psi_ind_rj_ik(rj, ik, deltak, deltaj, beta))
      } else{ # rj, ik, ij
        return(E_psi_ind_rj_ik_ij(rj, ik, ij, deltak, deltaj, beta))
      }
    }
    
    if(is.na(rj)){
      if(is.na(ik)){ # rk, ij
        return(E_psi_ind_rk_ij(rk, ij, deltak, deltaj, beta))
      } else{ # rk, ik, ij
        return(E_psi_ind_rk_ik_ij(rk, ik, ij, deltak, deltaj, beta))
      }
    }
    
    if(is.na(ij)){
      if(is.na(ik)){ # rk, rj
        return(E_psi_ind_rk_rj(rk, rj, deltak, deltaj, beta))
      } else{ # rk, rj, ik
        return(E_psi_ind_rk_rj_ik(rk, rj, ik, deltak, deltaj, beta))
      }
    }
    
    if(is.na(ik)){ # rk, rj, ij
      return(E_psi_ind_rk_rj_ij(rk, rj, ij, deltak, deltaj, beta))
    }
    
    return(E_psi_ind_rk_rj_ik_ij(rk, rj, ik, ij, deltak, deltaj, beta)) # rk, rj, ik, ij
  }
  
  # rk, rj
  E_psi_rk_rj <- function(rk, rj, deltak, deltaj, beta){
    if(rk > rj){
      out <- 1 - beta * deltaj / (deltaj + deltak) / (beta + deltak) * exp(- deltak * (rk - rj))
    } else{
      out <- deltak / (beta + deltak) + beta * deltak / (deltaj + deltak) / (beta + deltak) * exp(- deltaj * (rj - rk))
    }
    return(out)
  }
  
  E_psi_ind_rk_rj <- function(rk, rj, deltak, deltaj, beta){
    if(rk > rj){
      out <- deltaj * deltak / (deltaj + deltak) / (beta + deltak) * exp(- deltak * (rk - rj)) 
    } else{
      out <- deltaj * deltak / (deltaj + deltak) / (beta + deltak) * exp(- deltaj * (rj - rk)) 
    }
    return(out)
  }
  
  # ik, ij
  E_psi_ik_ij <- function(ik, ij, deltak, deltaj, beta){
    if(ik > ij){
      out <- 1
    } else{
      out <- exp( - (beta + deltak) * (ij - ik)) +
        deltak / (deltak + beta) * (1 - exp(- deltak * (ij - ik)))
    }
    return(out)
  }
  
  E_psi_ind_ik_ij <- function(ik, ij, deltak, deltaj, beta){
    if(ik > ij){
      out <- 0
    } else{
      out <- exp( - (beta + deltak) * (ij - ik))
    }
    return(out)
  }
  
  # rk, ij
  E_psi_rk_ij <- function(rk, ij, deltak, deltaj, beta){
    if(rk > ij){
      out <- deltak / (deltak + beta) * exp(- deltak * (rk - ij)) +
        (1 - exp(- deltak * (rk - ij)))
    } else{
      out <- deltak / (deltak + beta)
    }
    return(out)
  }
  
  E_psi_ind_rk_ij <- function(rk, ij, deltak, deltaj, beta){
    if(rk > ij){
      out <- deltak / (deltak + beta) * exp(- deltak * (rk - ij))
    } else{
      out <- 0
    }
    return(out)
  }
  
  # rj, ik
  E_psi_rj_ik <- function(rj, ik, deltak, deltaj, beta){
    if(rj > ik){
      out <- exp(- deltaj * (rj - ik)) +
        deltak / (deltak + beta) * (p2exp(rj - ik, deltak, deltaj) +
                                      pexp(rj - ik, deltak, lower.tail = F) +
                                      p2exp_cond(rj - ik, deltak, deltaj))
    } else{
      out <- 1
    }
    return(out)
  }
  
  E_psi_ind_rj_ik <- function(rj, ik, deltak, deltaj, beta){
    if(rj > ik){
      out <- deltak / (deltak + beta) * (pexp(rj - ik, deltak, lower.tail = F) +
                                           p2exp_cond(rj - ik, deltak, deltaj))
    } else{
      out <- 0
    }
    return(out)
  }
  
  # address
  
  # rk, rj, ij
  E_psi_rk_rj_ij <- function(rk, rj, ij, deltak, deltaj, beta){
    if(rk > ij){
      out <- deltak / (deltak + beta) * exp(- deltak * (rk - ij)) +
        (1 - exp(- deltak * (rk - ij)))
    } else{
      out <- deltak / (deltak + beta)
    }
    return(out)
    # if(rk > ij){
    #   if(rk > rj){
    #     out <- deltak / (deltak + deltaj) * exp(- deltak * (rk - rj)) +
    #       deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltak * (rk - rj))
    #     print('problem?')
    #     print(deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltak * (rk - rj)))
    #     print(deltak / (deltak + deltaj) * exp(- deltak * (rk - rj)))
    #   } else{
    #     out <- deltak / (deltak + deltaj) * exp(- deltaj * (rj - rk)) +
    #       deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltaj * (rj - rk))
    #     print('problem?')
    #     print(deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltaj * (rj - rk)))
    #     print(deltak / (deltak + deltaj) * exp(- deltaj * (rj - rk)))
    #   }
    # } else{
    #   #print('problem?')
    #   #print(deltak / (deltak + beta))
    #   out <- deltak / (deltak + beta)
    # }
    # return(out)
  }
  
  E_psi_ind_rk_rj_ij <- function(rk, rj, ij, deltak, deltaj, beta){
    if(rk > ij){
      out <- deltak / (deltak + beta) * exp(- deltak * (rk - ij))
    } else{
      out <- 0
    }
    return(out)
    # if(rk > ij){
    #   if(rk > rj){
    #     out <- deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltak * (rk - rj))
    #   } else{
    #     out <- deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltaj * (rj - rk))
    #   }
    # } else{
    #   out <- 0
    # }
    # return(out)
  }
  
  # address
  
  # rk, rj, ik
  E_psi_rk_rj_ik <- function(rk, rj, ik, deltak, deltaj, beta){
    if(rj < ik){
      u <- 1
    } else{
      u <- exp(- deltaj * (rj - ik))
    }
    if(rj < rk){
      v <- 0
    } else{
      v <- exp(- beta * (rk - ik)) * (1 - exp(- deltaj * (rj - rk)))
    }
    if(rj < ik){
      w <- 0
    } else if(rk < rj){
      w <- deltak / (deltak + beta) * (1 - exp(- deltaj * (rk - ik))) * exp(- deltaj * (rj - rk))
    } else{
      w <- deltak / (deltak + beta) * (1 - exp(- deltaj * (rj - ik)))
    }
    return(u+v+w)
    
    # if(ik > rj){
    #   out <- 1
    # } else{
    #   if(rk > rj){
    #     out <- deltak / (deltak + deltaj) * exp(- deltak * (rk - rj)) +
    #       deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltak * (rk - rj))
    #   } else{
    #     out <- deltak / (deltak + deltaj) * exp(- deltaj * (rj - rk)) + 
    #       deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltaj * (rj - rk)) +
    #       exp(- beta * (rk - ik)) * (1 - exp(- deltaj * (rj - rk)))
    #   }
    # }
    # return(out)
  }
  
  E_psi_ind_rk_rj_ik <- function(rk, rj, ik, deltak, deltaj, beta){
    if(rj < ik){
      w <- 0
    } else if(rk < rj){
      w <- deltak / (deltak + beta) * (1 - exp(- deltaj * (rk - ik))) * exp(- deltaj * (rj - rk))
    } else{
      w <- deltak / (deltak + beta) * (1 - exp(- deltaj * (rj - ik)))
    }
    return(w)
    # if(ik > rj){
    #   out <- 0
    # } else{
    #   if(rk > rj){
    #     out <- deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltak * (rk - rj))
    #   } else{
    #     out <- deltak / (deltak + beta) * deltaj / (deltaj + deltak) * exp(- deltaj * (rj - rk))
    #   }
    # }
    # return(out)
  }
  
  # rk, ik, ij
  E_psi_rk_ik_ij <- function(rk, ik, ij, deltak, deltaj, beta){
    taukj <- min(rk, ij) - min(ik, ij)
    out <- exp(- beta * taukj)
    return(out)
  }
  
  E_psi_ind_rk_ik_ij <- function(rk, ik, ij, deltak, deltaj, beta){
    out <- 0
    if(ik < ij){
      if(ij < rk){
        out <- exp( - beta * (ij - ik))
      }
    }
    return(out)
  }
  
  # rj, ik, ij
  E_psi_rj_ik_ij <- function(rj, ik, ij, deltak, deltaj, beta){
    if(ik > ij){
      out <- 1
    } else{
      out <- deltak / (deltak + beta) * (1 - exp(- deltak * (ij - ik))) +
        exp(- (beta + deltak) * (ij - ik))
    }
    return(out)
  }
  
  E_psi_ind_rj_ik_ij <- function(rj, ik, ij, deltak, deltaj, beta){
    if(ik > ij){
      out <- 0
    } else{
      out <- exp(- (beta + deltak) * (ij - ik))
    }
    return(out)
  }
  
  # rk, rj, ik, ij
  E_psi_rk_rj_ik_ij <- function(rk, rj, ik, ij, deltak, deltaj, beta){
    taukj <- min(rk, ij) - min(ik, ij)
    out <- exp(- beta * taukj)
    return(out)
  }
  
  E_psi_ind_rk_rj_ik_ij <- function(rk, rj, ik, ij, deltak, deltaj, beta){
    out <- 0
    if(ik < ij){
      if(ij < rk){
        out <- exp( - beta * (ij - ik))
      }
    }
    return(out)
  }
  
  p2exp <- function(z, lambda2, lambda1){
    if(lambda2 != lambda1){ # hypoexponential
      return((lambda1 * exp(- lambda2 * z) - lambda2 * exp(- lambda1 * z) + lambda2 - lambda1) / 
               (lambda2 - lambda1))
    } else{ # gamma
      return(pgamma(z, shape = 2, rate = lambda2))
    }
  }
  
  p2exp_cond <- function(z, lambda2, lambda1){
    lambda <- lambda2 - lambda1
    if(lambda != 0){ # hypoexponential
      out <- exp(- lambda1 * z) * 
        (lambda2 / lambda * (1 - exp(- lambda * z)) - (1 - exp(- lambda2 * z)))
    } else{ # not hypoexponential
      out <- exp(- lambda1 * z) * (lambda1 * z + 1 - exp(- lambda1 * z))
    }
    return(out)
  }
  
  ##### main function #####
  
  # initialize
  n <- length(r)
  beta <- beta / N
  r1 <- min(i, r, na.rm = T)
  
  # change of variable to delta
  if(n < N){
    B <- beta * (N - n)
    delta <- gamma + B
  } else{ # handles entire population infected
    if(n == N){delta = gamma}
  }
  
  # calculate log likelihood (line six)
  ind <- order(r)
  ind <- ind[is.na(i)]
  ind <- ind[1:min(A,length(ind))]
  ip = - delta * (r[ind] - r1)
  ia = rep(-log(length(ind)), length(ind))
  z = ia + ip
  
  # copy and paste from is.integer documentation
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }
  
  if((any(beta < 0)) | (any(gamma < 0)) |
     (!is.wholenumber(N)) | (N <= 0) |
     (!is.wholenumber(A)) | (A <= 0)){
    # invalid parameters
    return(1e15)
  } else{
    
    # evaluate psi and chi terms
    psichi <- rep(0, n)
    for(j in (1:n)){
      X <- 0
      Y <- 0
      rj <- r[j]
      ij <- i[j]
      for(k in (1:n)[-j]){
        rk <- r[k]
        ik <- i[k]
        x <- E_psi_ind(rk, rj, ik, ij, delta, delta, beta)
        y <- E_psi(rk, rj, ik, ij, delta, delta, beta)
        #print(c(k,j,log(y)))
        X = X + beta * x / y
        Y = Y + log(y)
      }
      #print(c(j, log(X), Y))
      psichi[j] = Y + log(X)
    }
    
    # line eight
    z <- as.vector(z)
    for(alpha in 1:length(ind)){z[alpha] <- z[alpha] + sum(psichi[-ind[alpha]])}
    z <- matrixStats::logSumExp(z)
    a <- n * log(gamma / delta)
    
    # negative log likelihoods
    return(-(a + z))  
  }
}


