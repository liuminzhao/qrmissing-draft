## Time-stamp: <liuminzhao 03/29/2013 01:03:23>
## WRAP UP BiMLESigma.f

dyn.load('BiMLESigma.so')
BiQRGradient <- function(y, R, x, tau, niter = 1000){
  n <- length(R)
  param <- rep(0, 14)
  paramsave <- matrix(0, niter, 15)
  mod <- .Fortran("BiQRGradientf",
                  y = as.double(y),
                  R = as.integer(R),
                  x = as.double(x),
                  tau = as.double(tau),
                  n = as.integer(n),
                  niter = as.integer(niter),
                  param = as.double(param),
                  paramsave = as.double(paramsave)
                  )
  return(mod)
}
