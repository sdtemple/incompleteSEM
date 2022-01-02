# incompleteSEM Simulation Study
# Seth Temple, sdtemple@uw.edu

# Setup -------------------------------------------------------------------

library(incompleteSEM)

# epidemic model parameters
N <- 300
beta <- 2.5
gamma <- 1
rates <- c(beta, gamma)
init <- c(2,2)

# simulation parameters
simCt <- 0
simSize <- 500
simMin <- 30
simP <- c(1/5,2/5,3/5,4/5)
simQ <- 1
simPrint <- 1

# storage parameters
simStore <- array(NA, dim = c(simSize,2,2+length(simP)))
simI <- rep(NA, simSize)

# Simulation Study --------------------------------------------------------

while(simCt < simSize){
  # simulate epidemic
  epi <- rgsem(rates, N)
  epi <- filter_gsem(epi)
  while(dim(epi)[1] < simMin){
    epi <- rgsem(rates, N)
    epi <- filter_gsem(epi)
  }
  simCt <- simCt + 1
  
  # record epidemic size
  simI[simCt] <- dim(epi)[1]

  # complete epidemic
  r <- epi[,2]
  i <- epi[,1]
  simStore[simCt,1:2,1] <- mle_complete_gsem(r, i, N)
  #simStore[simCt,3,1] <- simStore[simCt,1,1] / simStore[simCt,2,1]
  
  # partial epidemic
  epi <- sort_gsem(epi)
  r <- epi[,2]
  simStore[simCt,1:2,2] <- nlm(pbla_partial_gsem, init, r, N)$estimate
  #simStore[simCt,3,2] <- simStore[simCt,1,2] / simStore[simCt,2,2]
  
  # incomplete epidemics
  for(j in 1:length(simP)){
    p <- simP[j]
    iepi <- decomplete_gsem(epi, p, simQ)
    iepi <- sort_gsem(iepi)
    r <- iepi[,2]
    i <- iepi[,1]
    simStore[simCt,1,2+j] <- nlm(pbla_incomplete_beta, init[1], r, i, N)$estimate
    simStore[simCt,2,2+j] <- mle_gamma(r, i)[1]
    # simStore[simCt,1:2,2+j] <- nlm(pbla_incomplete_gsem, init, r, i, N)$estimate
    #simStore[simCt,3,2+j] <- simStore[simCt,1,2+j] / simStore[simCt,2,2+j]
  }
  
  # printing
  if(!(simCt %% simPrint)){print(simCt)}
}

# saving
fn <- paste('../data/ii-', beta, '-', gamma, '-', N, '-mat.rds', sep = '')
saveRDS(simStore, fn)
fn <- paste('../data/ii-', beta, '-', gamma, '-', N, '-vec.rds', sep = '')
saveRDS(simI, fn)
