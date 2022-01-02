#' MLE for complete general stochastic epidemic model
#'
#' Compute MLE for complete epidemic observations.
#'
#' @param r numeric vector: removal times
#' @param i numeric vector: infection times
#' @param N integer: population size
#'
#' @return MLE for (beta, gamma)
#'
#' @export
mle_complete_gsem = function(r, i, N){
    n = length(r)
    t = 0
    for(j in 1:n){
        t = t + sum(sapply(i, min, r[j]) - sapply(i, min, i[j]))
    }
    ri = sum(r - i)
    g = n / ri
    b = (n - 1) / (t + (N - n) * ri) * N
    return(c(b,g))
}