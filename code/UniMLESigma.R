## Time-stamp: <liuminzhao 03/28/2013 10:50:54>
## WRAP UP UniMLESigma.f

rm(list = ls())
dyn.load('UniMLESigma.so')
QRGradient <- function(y, S, x, tau, niter = 1000){
  n <- length(y)
  param <- rep(0, 7)
  paramsave <- matrix(0, niter, 8)
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
  a <- matrix(mod$paramsave, mod$niter, 8)
  colnames(a) <- c('Gamma0', 'Gamma1', 'Beta0', 'Beta1', 'Sigma1', 'Sigma0', 'phi', 'NLL')
  for (i in 1:8){
    plot(ts(a[, i]), main = colnames(a)[i])
  }
}
