#' @title Estimate expected homogeneity index, or required values of mutation rate or population for a user defined homogeneity index.
#' @params mu Vector of mutation value. Default is NULL.
#' @params N Vector of population sizes. Default is NULL.
#' @params H Vector of homogeneity index. Default is NULL.
#' @details Computes value of H/N/mu based on user definedc
#' @return A data.frame with given an expected (under equilibrium) estimates of mutation rate (mu), population size (N), and homogeneity (H).
#' @examples 
#' mu=c(0.01,0.01)
#' H=c(0.8,0.5)
#' muNH(mu=mu,H=H)
#' N=c(1000,1000)
#' H=c(0.8,0.5)
#' muNH(H=H,N=N)
#' N=c(1000,1000)
#' muc(0.01,0.001)
#' muNH(mu=mu,N=N)
muNH <- function(mu=NULL,N=NULL,H=NULL)
{
  if (any(is.null(mu))&!is.null(N)&!is.null(H))
  {
    mu = (-H + 1)/(2*H*N)
  }
  if (any(is.null(N))&!is.null(H)&!is.null(H))
  {
    N = (-H + 1)/(2*H*mu)
  }
  if (any(is.null(H))&!is.null(mu)&!is.null(N))
  {
    H = 1/(2*N*mu + 1)
  }
  return(data.frame(N,mu,H))
}