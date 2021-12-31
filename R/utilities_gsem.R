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

#' Missing epidemic data
#' 
#' Impute NAs for infection and removal times
#' 
#' @param epi infection and removal times
#' @param p expected proportion of complete pairs observed
#' @param q probability infection time missing
impute_gsem <- function(epi, p, q = 1){
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