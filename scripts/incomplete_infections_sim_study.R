# incompleteSEM Simulation Study
# Seth Temple, sdtemple@uw.edu

# Setup -------------------------------------------------------------------

library(incompleteSEM)

# epidemic model parameters
N <- 400
beta <- 1.8
gamma <- 1
rates <- c(beta, gamma)

# simulation parameters
simCt <- 0
simSize <- 200
simMin <- 10
simP <- c(1/4,1/2,3/4)
simQ <- 1
simPrint <- 10

# storage parameters
simStore <- array(NA, dim = c(simSize,3,2+length(simP)))
simI <- rep(NA, simSize)

# testing -----------------------------------------------------------------

# epi <- rgsem(c(beta, gamma), N)
# epi <- filter_gsem(epi)
# 
# likelihood_complete_gsem(rates, epi[,2], epi[,1], N)
# mle_complete_gsem(epi[,2], epi[,1], N)
# pbla_incomplete_gsem(rates, epi[,2], epi[,1], N)
# 
# pbla_partial_gsem(rates, epi[,2], N)
# 
# epi <- sort_gsem(epi)
# 
# epi2 <- impute_gsem(epi, 1/2)
# pbla_incomplete_gsem(rates, epi2[,2], epi2[,1], N)
# 
# 
# nlm(pbla_incomplete_gsem, c(1,1), r=epi2[,2], i=epi2[,1], N=N)$estimate
# nlm(pbla_partial_gsem, c(1,1), r=epi2[,2], N=N)$estimate
# mle_complete_gsem(epi[,2], epi[,1], N)

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
  simStore[simCt,3,1] <- simStore[simCt,1,1] / simStore[simCt,2,1]
  
  # partial epidemic
  epi <- sort_gsem(epi)
  r <- epi[,2]
  simStore[simCt,1:2,2] <- nlm(pbla_partial_gsem, rates, r, N)$estimate
  simStore[simCt,3,2] <- simStore[simCt,1,2] / simStore[simCt,2,2]
  
  # incomplete epidemics
  for(j in 1:length(simP)){
    p <- simP[j]
    iepi <- decomplete_gsem(epi, p, simQ)
    iepi <- sort_gsem(iepi)
    r <- iepi[,2]
    i <- iepi[,1]
    simStore[simCt,1:2,2+j] <- nlm(pbla_incomplete_gsem, rates, r, i, N)$estimate
    simStore[simCt,3,2+j] <- simStore[simCt,1,2+j] / simStore[simCt,2,2+j]
  }
  
  # printing
  if(!(simCt %% simPrint)){print(simCt)}
}

# saving
fn <- paste('../data/ii-', beta, '-', gamma, '-', N, '-mat.rds', sep = '')
saveRDS(simStore, fn)
fn <- paste('../data/ii-', beta, '-', gamma, '-', N, '-vec.rds', sep = '')
saveRDS(simI, vec)
