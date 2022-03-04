#' @title Function that computes sumlogs
#'
#' @param m a vector
#'
#' @return the result of sumlogs

#'
#' @examples \dontrun{
#' library(rmSCGLR)
#'
#' sumlogs(c(1,2,3))
#'
#' }
sumlogs <- function(m){
  M <- max(m)
  return(M+log(sum(exp(m-M))))
}

wtScale <-function(x,w) {
  xc <- x-sum(w*x)
  v <- sum(xc^2*w)
  xcr <- xc/sqrt(v)
  return(xcr)
}
