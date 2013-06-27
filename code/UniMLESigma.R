##' Time-stamp: <liuminzhao 06/27/2013 14:01:32>
##' WRAP UP UniMLESigma.f
##' 2013/06/26 heter2 as default method
##' dyn load

## dyn.load('UniMLESigma.so')
## dyn.load('UniMLESigmaHeter1.so')
## dyn.load('UniMLESigmaHeter2.so')
QRGradient <- function(y, S, x, tau, niter = 1000, method = "heter2"){
  n <- length(y)
  if (method == "homo"){
    param <- rep(0, 7)
    paramsave <- matrix(0, niter, 8)
    mod <- .Fortran("qrgradientf",
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
    mod <- .Fortran("qrgradientheter1f",
                    y = as.double(y),
                    S = as.integer(S),
                    x = as.double(x),
                    tau = as.double(tau),
                    n = as.integer(n),
                    niter = as.integer(niter),
                    param = as.double(param),
                    paramsave = as.double(paramsave)
                    )
  } else if (method == "heter2") {
    if (!is.loaded('qrgradienth2f')) {
      if (file.exists('UniMLESigmaHeter2.so')) {
        dyn.load('UniMLESigmaHeter2.so')
      } else {
        if (file.exists('~/Documents/qrmissing/code/UniMLESigmaHeter2.so')) {
          dyn.load("~/Documents/qrmissing/code/UniMLESigmaHeter2.so")
        } else stop('no shared library found.')
      }
    }

    param <- rep(0, 9)
    paramsave <- matrix(0, niter, 10)
    mod <- .Fortran("QRGradientH2f",
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
  class(mod) <- "QRGradient"
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
  } else if (mod$method == "heter2") {
    a <- matrix(mod$paramsave, mod$niter, 10)
    colnames(a) <- c('Gamma0', 'Gamma1', 'Beta0', 'Beta1', 'alpha01', 'ALpha11', 'alpha02', 'alpha12','phi', 'NLL')
    for (i in 1:10){
      plot(ts(a[, i]), main = colnames(a)[i])
    }
  }
}

print.QRGradient <- function(mod, ...) {
  cat('Coefficients:\n')
  print(coef(mod))
}

summary.QRGradient <- function(mod, ...){
  n <- mod$n
  S <- mod$S
  tau <- mod$tau
  param <- mod$param

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(S)/n, '\n')
  cat('Estimated pi:', param[9], '\n')
  cat('Quantile: ', tau, '\n')
  cat('Quantile regression coefficients: \n')
  mycoef <- matrix(param[1:8], 4, 2, byrow = TRUE)
  colnames(mycoef) <- c('Intercept', 'Slope')
  rownames(mycoef) <- c('Gamma', 'Beta', 'Alpha1', 'Alpha0')
  print(mycoef)
}

plot.QRGradient <- function(mod, ...){
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
  } else if (mod$method == "heter2") {
    a <- matrix(mod$paramsave, mod$niter, 10)
    colnames(a) <- c('Gamma0', 'Gamma1', 'Beta0', 'Beta1', 'alpha01', 'ALpha11', 'alpha02', 'alpha12','phi', 'NLL')
    for (i in 1:10){
      plot(ts(a[, i]), main = colnames(a)[i])
    }
  }
}

coef.QRGradient <- function(mod, ...){
  ans <- matrix(mod$param[1:2], 1, 2)
  colnames(ans) <- c('Intercept', 'Slope')
  return(ans)
}
