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