## Time-stamp: <liuminzhao 03/27/2013 22:29:52>
## WRAP UP UniMLESigma.f

rm(list = ls())
dyn.load('UniMLESigma.so')
QRGradient <- function(y, S, x, tau, niter = 1000){
  n <- length(y)
  param <- rep(0, 7)
  paramsave <- matrix(0, 8, niter)
  mod <- .Fortran("QRGradientf",
                  y = as.double(y),
                  S = as.integer(S),
                  x = as.double(x),
                  tau = as.double(tau),
                  n = as.integer(n),
                  niter = as.integer(niter),
                  param = as.double(param),
                  paramsave = as.double(paramsave)
                  )
  return(mod)
}

Diagnose <- function(mod){
  a <- matrix(mod$paramsave, 8, mod$niter)
  for (i in 1:8){
    plot(ts(a[i,]), main = rownames(a)[i])
  }
}
