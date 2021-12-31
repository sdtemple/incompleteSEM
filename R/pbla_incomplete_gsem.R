#' PBLA for incomplete general stochastic epidemic model
#' 
#' Compute approximate likelihood for incomplete epidemic observations.
#' 
#' @param rates numeric vector: rates (beta, gamma)
#' @param r numeric vector: removal times
#' @param i numeric vector: of infection times
#' @param N integer: population size
#' @param A integer: plausible patient zeros
#' 
#' @return negative log likelihood
#' 
#' @export
pbla_incomplete_gsem <- function(rates, r, i, N, A = 1){
  
  beta <- rates[1]
  gamma <- rates[2]
  
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
    #print(out)
    return(out)
  }
  
  E_psi_ind_rk_rj <- function(rk, rj, deltak, deltaj, beta){
    if(rk > rj){
      out <- deltaj * deltak / (deltaj + deltak) / (beta + deltak) * exp(- deltak * (rk - rj)) 
    } else{
      out <- deltaj * deltak / (deltaj + deltak) / (beta + deltak) * exp(- deltaj * (rj - rk)) 
    }
    #print(out)
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
    #print(out)
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
    #print(out)
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
    #print(out)
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
  
  # unpleasant integral
  # unpleasant_integral <- function(b, a, rk, rj, ik, deltak, deltaj, beta){
  #   outside <- - deltaj * exp(beta * ik) * exp(- deltaj * rj) / (beta - deltaj)
  #   inside <- exp(- (beta - deltaj) * b) - exp(- (beta - deltaj) * a)
  #   return(inside * outside)
  # }
  
  # unpleasant integral
  unpleasant_integral <- function(rk, rj, ik, deltak, deltaj, beta){
    a <- max(rj - rk, 0)
    b <- rj - ik
    propto <- 1 / (pexp(b) - pexp(a))
    outside <- exp(- beta * (rj - ik))
    inside <- deltaj / (deltaj - beta) * (exp(- (deltaj - beta) * a) - exp(- (deltaj - beta) * b))
    return(inside * outside * propto)

    # outside <- - beta * (rj - ik)
    # inside <- log(deltaj) - log(deltaj - beta) + log(exp(- (deltaj - beta) * a) - exp(- (deltaj - beta) * b))
    # propto <- 1 / (pexp(b) - pexp(a))
    # return(inside + outside + log(propto))
    
    #return(log(deltaj) - log(deltaj + beta))
  }
  
  # rk, rj, ik
  E_psi_rk_rj_ik <- function(rk, rj, ik, deltak, deltaj, beta){
    #print(deltaj)
    if(rj < ik){
      u <- 1
      #u <- log(1)
    } else{
      u <- exp(- deltaj * (rj - ik))
      #u <- - deltaj * (rj - rk)
    }
    if(rj < rk){
      v <- 0
    } else{
      v <- exp(- beta * (rk - ik)) * (1 - exp(- deltaj * (rj - rk)))
      # v <- - beta * (rk - ik) + log(1 - exp(- deltaj * (rj - rk)))
      # v <- exp(v)
    }
    if(rj < ik){
      w <- 0
    } else if(rk < rj){
      #w <- unpleasant_integral(rk, ik, rk, rj, ik, deltak, deltaj, beta) *
      w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) *
       (1 - exp(- deltaj * (rk - ik))) * exp(- deltaj * (rj - rk))
      # w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) + 
      #   log(1 - exp(- deltaj * (rk - ik))) - 
      #   deltaj * (rj - rk)
      # w <- exp(w)
    } else{
      #w <- unpleasant_integral(rk, ik, rk, rj, ik, deltak, deltaj, beta) *
      w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) *
       (1 - exp(- deltaj * (rj - ik)))
      # w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) + 
      #   log(1 - exp(- deltaj * (rj - ik)))
      # w <- exp(w)
    }
    #print(u+v+w)
    #print('u')
    #print(u)
    #print('v')
    #print(v)
    #if(v < 0.9){print(v)}
    #if(v == 0){print(c(u,w))}
    #print('w')
    #print(w)
    #print(u+v+w)
    #print(u)
    #print(v)
    #print(w)
    return(u+v+w)
    #return(u+v)
    
    # if(rj < ik){
    #   u <- 1
    # } else{
    #   u <- exp(- deltaj * (rj - ik))
    # }
    # if(rj < rk){
    #   v <- 0
    # } else{
    #   v <- exp(- beta * (rk - ik)) * (1 - exp(- deltaj * (rj - rk)))
    # }
    # if(rj < ik){
    #   w <- 0
    # } else if(rk < rj){
    #   w <- deltak / (deltak + beta) * (1 - exp(- deltaj * (rk - ik))) * exp(- deltaj * (rj - rk))
    # } else{
    #   w <- deltak / (deltak + beta) * (1 - exp(- deltaj * (rj - ik)))
    # }
    # return(u+v+w)
    
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
      #w <- unpleasant_integral(rk, ik, rk, rj, ik, deltak, deltaj, beta) *
      w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) *
       (1 - exp(- deltaj * (rk - ik))) * exp(- deltaj * (rj - rk))
      # w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) + 
      #   log(1 - exp(- deltaj * (rk - ik))) - 
      #   deltaj * (rj - rk)
      # w <- exp(w)
    } else{
      #w <- unpleasant_integral(rk, ik, rk, rj, ik, deltak, deltaj, beta) *
      w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) *
       (1 - exp(- deltaj * (rj - ik)))
      # w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) + 
      #   log(1 - exp(- deltaj * (rj - ik)))
      # w <- exp(w)
    }
    #   w <- 0
    # } else if(rk < rj){
    #   #w <- unpleasant_integral(rk, ik, rk, rj, ik, deltak, deltaj, beta) *
    #   w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) * 
    #     (1 - exp(- deltaj * (rk - ik))) * exp(- deltaj * (rj - rk)) 
    # } else{
    #   #w <- unpleasant_integral(rk, ik, rk, rj, ik, deltak, deltaj, beta) *
    #   w <- unpleasant_integral(rk, rj, ik, deltak, deltaj, beta) * 
    #     (1 - exp(- deltaj * (rj - ik)))
    # }
    #print(w)
    return(w)
    
    # if(rj < ik){
    #   w <- 0
    # } else if(rk < rj){
    #   w <- deltak / (deltak + beta) * (1 - exp(- deltaj * (rk - ik))) * exp(- deltaj * (rj - rk))
    # } else{
    #   w <- deltak / (deltak + beta) * (1 - exp(- deltaj * (rj - ik)))
    # }
    # return(w)
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
  
  ##### main function ######
  
  # case for complete likelihood
  if(!any(is.na(i))){
    if(!any(is.na(r))){
      return(likelihood_complete_gsem(r, i, beta, gamma, N))
    }
  }
  
  # case for partially observed likelihood
  if(all(is.na(i))){
    return(pbla_partial_gsem(r, beta, gamma, N))
  }
  
  # set up
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
  if(any(is.na(i))){
    ind <- order(r)
    ind <- ind[is.na(i)]
    ind <- ind[1:min(A,length(ind))]
    ip <- - delta * (r[ind] - r1)
    ia <- rep(-log(length(ind)), length(ind))
    z <- ia + ip
  } else{
    ind <- which.min(i)
    z <- 0
  }
  
  # copy and paste from is.integer documentation
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }
  
  if((any(beta < 0)) | (any(gamma < 0)) |
     (!is.wholenumber(N)) | (N <= 0) |
     (!is.wholenumber(length(ind))) | (length(ind) <= 0)){
    # invalid parameters
    return(1e15)
  }
  
  # index subsets
  iincomplete <- which(is.na(i), arr.ind = T)
  rincomplete <- which(is.na(r), arr.ind = T)
  complete <- which(!is.na(r) & !is.na(i), arr.ind = T)
  incomplete <- c(iincomplete, rincomplete)
  
  #print(iincomplete)
  #print(rincomplete)
  #print(complete)
  #print(incomplete)
  
  ### outside of the integral ###
  
  log_phi_complete <- - beta * (N - n) * sum(r[complete] - i[complete]) 
  #print(log_phi_complete)
  
  log_f_complete <- - gamma * sum(r[complete] - i[complete]) + length(complete) * log(gamma)
  #print(log_f_complete)
  
  log_psi_complete <- rep(0, n)
  for(k in 1:length(complete)){
    rk <- r[complete[k]]
    ik <- i[complete[k]]
    for(j in (1:length(complete))[-k]){
      ij <- i[complete[j]]
      tau <- min(rk, ij) - min(ik, ij)
      log_psi_complete[complete[j]] <- log_psi_complete[complete[j]] + tau
    }
  }
  log_psi_complete <- - beta * log_psi_complete
  #print(log_psi_complete)
  
  ### inside of the integral ###
  
  psichi <- rep(0, n)
  if(length(complete) > 0){
    # completes and incompletes
    for(j in 1:length(complete)){
      X <- 0
      Y <- 0
      rj <- r[complete[j]]
      ij <- i[complete[j]]
      for(k in 1:length(incomplete)){
        rk <- r[incomplete[k]]
        ik <- i[incomplete[k]]
        y <- E_psi(rk, rj, ik, ij, delta, delta, beta)
        x <- E_psi_ind(rk, rj, ik, ij, delta, delta, beta)
        #print(c(incomplete[k], complete[j], x))
        #if(y < .9){print(c(incomplete[k], complete[j], y))}
        X <- X + beta * min(x / y, 1)
        y <- min(y, 1)
        Y <- Y + log(y)
      }
      #print(c(complete[j], Y, log(X)))
      for(k in (1:length(complete))[-j]){
        rk <- r[complete[k]]
        ik <- i[complete[k]]
        X <- X + as.numeric((ik < ij) & (ij < rk)) * beta
      }
      #print(c(complete[j], Y, log(X)))
      psichi[complete[j]] <- Y + log(X)
    }
    for(j in 1:length(incomplete)){
      X <- 0
      Y <- 0
      rj <- r[incomplete[j]]
      ij <- i[incomplete[j]]
      for(k in 1:length(complete)){
        rk <- r[complete[k]]
        ik <- i[complete[k]]
        y <- E_psi(rk, rj, ik, ij, delta, delta, beta)
        x <- E_psi_ind(rk, rj, ik, ij, delta, delta, beta)
        
        #if(y < .9){print(c(complete[k], incomplete[j], rk, rj, ik, ij, y))}
        #print(c(complete[k],incomplete[j],x))
        X <- X + beta * min(x / y, 1)
        y <- min(y, 1)
        #print(c(complete[k], incomplete[j], rk, rj, ik, ij, y))
        Y <- Y + log(y)
      }
      #print(c(incomplete[j], ij, Y, log(X)))
      for(k in (1:length(incomplete))[-j]){
        rk <- r[incomplete[k]]
        ik <- i[incomplete[k]]
        y <- E_psi(rk, rj, ik, ij, delta, delta, beta)
        x <- E_psi_ind(rk, rj, ik, ij, delta, delta, beta)
        X <- X + beta * min(x / y, 1)
        y <- min(y, 1)
        Y <- Y + log(y)
      }
      #print(c(incomplete[j], Y, log(X)))
      psichi[incomplete[j]] <- Y + log(X)
    }
  } else{
    # incompletes only
    for(j in 1:length(incomplete)){
      X <- 0
      Y <- 0
      rj <- r[incomplete[j]]
      ij <- i[incomplete[j]]
      for(k in (1:length(incomplete))[-j]){
        rk <- r[incomplete[k]]
        ik <- i[incomplete[k]]
        y <- E_psi(rk, rj, ik, ij, delta, delta, beta)
        x <- E_psi_ind(rk, rj, ik, ij, delta, delta, beta)
        X <- X + beta * min(x / y, 1)
        y <- min(y, 1)
        Y <- Y + log(y)
      }
      #print(c(Y, log(X)))
      psichi[incomplete[j]] <- Y + log(X)
    }
  }
  
  #print(psichi)
  #print(sum(psichi))
  
  # line eight
  z <- as.vector(z)
  for(alpha in 1:length(ind)){
    z[alpha] <- z[alpha] + 
      sum(psichi[-ind[alpha]]) +
      sum(log_psi_complete[-ind[alpha]])
  }
  z <- matrixStats::logSumExp(z)
  a <- length(incomplete) * log(gamma / delta)
  
  #print(-log_phi_complete)
  #print(-sum(log_psi_complete))
  #print(-log_f_complete)
  
  # negative log likelihoods
  return(-(a + z + log_phi_complete + log_f_complete))
}
