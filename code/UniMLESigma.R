## Time-stamp: <liuminzhao 04/13/2013 22:36:43>
## WRAP UP UniMLESigma.f
dyn.load('UniMLESigma.so')
dyn.load('UniMLESigmaHeter1.so')
QRGradient <- function(y, S, x, tau, niter = 1000, method = "homo"){
  n <- length(y)
  if (method == "homo"){
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
  } else if (method == "heter1") {
    param <- rep(0, 9)
    paramsave <- matrix(0, niter, 10)
    mod <- .Fortran("QRGradientHeter1f",
                    y = as.double(y),
                    S = as.integer(S),
                    x = as.double(x),
                    tau = as.double(tau),
                    n = as.integer(n),
                    niter = as.integer(niter),
                    param = as.double(param),
                    paramsave = as.double(paramsave)
                    )
  }

  mod$method <- method
  return(mod)
}

Diagnose <- function(mod){
  if (mod$method == "homo") {
    a <- matrix(mod$paramsave, mod$niter, 8)
    colnames(a) <- c('Gamma0', 'Gamma1', 'Beta0', 'Beta1', 'Sigma1', 'Sigma0', 'phi', 'NLL')
    for (i in 1:8){
      plot(ts(a[, i]), main = colnames(a)[i])
    }
  } else if (mod$method == "heter1") {
    a <- matrix(mod$paramsave, mod$niter, 10)
    colnames(a) <- c('Gamma0', 'Gamma1', 'Beta0', 'Beta1', 'Sigma1', 'Sigma0', 'phi', 'NLL','Heter1', 'Heter2')
    for (i in 1:10){
      plot(ts(a[, i]), main = colnames(a)[i])
    }
  }
}
