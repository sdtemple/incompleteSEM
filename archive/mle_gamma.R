#' MLE for Removal Rate
#' 
#' Compute the maximum likelihood estimate for the removal rate.
#' 
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' 
#' @return named vector (estimate, sample size)
#'  
#' @export 
mle_gamma <- function(r, i){
  ind <- (!is.na(r)) * (!is.na(i))
  ind <- which(ind == 1, arr.ind=T)
  r <- r[ind]
  i <- i[ind]
  ri <- r - i
  return(c(mle = length(ri) / sum(ri), n = length(ri)))
}
