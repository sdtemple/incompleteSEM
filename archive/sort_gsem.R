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