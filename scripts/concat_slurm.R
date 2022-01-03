# Concat slurm array outputs
# Seth Temple, sdtemple@uw.edu

flavor <- 'peanut'
beta <- 4
gamma <- 1
N <- 300
indices <- 1:40

vec <- c()
for(i in indices){
  folder <- paste(flavor, beta, gamma, N, sep = '-')
  file <- paste(flavor, beta, gamma, N, i, 'vec', sep = '-')
  file <- paste(file, '.rds', sep = '')
  file <- paste(folder, file, sep = '/')
  vec <- c(vec, readRDS(file))
}
saveRDS(vec, paste(folder, '-vec.rds', sep = ''))

library(abind)
mat <- array(NA, dim = c(0,2,6))
for(i in indices){
  folder <- paste(flavor, beta, gamma, N, sep = '-')
  file <- paste(flavor, beta, gamma, N, i, 'mat', sep = '-')
  file <- paste(file, '.rds', sep = '')
  file <- paste(folder, file, sep = '/')
  mat <- abind(mat, readRDS(file), along=1)
}
saveRDS(mat, paste(folder, '-mat.rds', sep = ''))
