## Time-stamp: <liuminzhao 03/27/2013 19:58:10>
## WRAP UP UniMLESigma.f

rm(list = ls())
dyn.load('UniMLESigma.so')
QRGradient <- function(y, S, x, tau, niter = 400){
  n <- length(y)
  param <- rep(0, 7)
  mod <- .Fortran("QRGradientf",
                  y = as.double(y),
                  S = as.integer(S),
                  x = as.double(x),
                  tau = as.double(tau),
                  n = as.integer(n),
                  niter = as.integer(niter),
                  param = as.double(param)
                  )
  return(mod)
}
