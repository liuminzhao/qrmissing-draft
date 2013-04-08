## Time-stamp: <liuminzhao 04/04/2013 18:32:32>
## WRAP UP UniMLESigma.f

dyn.load('UniMLESigma.so')
QRGradient <- function(y, S, x, tau, niter = 1000){
  n <- length(y)
  param <- rep(0, 9)
  paramsave <- matrix(0, niter, 10)
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
  a <- matrix(mod$paramsave, mod$niter, 10)
  colnames(a) <- c('Gamma0', 'Gamma1', 'Beta0', 'Beta1', 'Sigma1', 'Sigma0', 'phi', 'NLL',
                   'Heter1', 'Heter2')
  for (i in 1:10){
    plot(ts(a[, i]), main = colnames(a)[i])
  }
}
