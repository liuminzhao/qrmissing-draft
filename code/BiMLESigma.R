## Time-stamp: <liuminzhao 07/26/2013 12:50:54>
## WRAP UP BiMLESigma.f

BiQRGradient <- function(y, R, X, tau=0.5, niter = 1000, sp = rep(0, 1 + 2*dim(X)[2]),
                         method = 'heter2', init = NULL, tol = 0.00001){
  ## LOAD SHARED LIBRARY
  if (method == 'heter2') {
    if (!is.loaded('biqrgradienth2f')) {
      if (file.exists('BiMLESigmaH2.so')) {
        dyn.load('BiMLESigmaH2.so')
      } else {
        if (file.exists('~/Documents/qrmissing/code/BiMLESigmaH2.so')) {
          dyn.load("~/Documents/qrmissing/code/BiMLESigmaH2.so")
        } else stop('no shared library found.')
      }
    }
  }

  n <- length(R)
  if (method == 'homo') {
    param <- rep(0, 14)
    param[1] = 0
    param[2] = 0
    param[3] = 0
    param[4] = 0
    param[5] = 1
    param[6] = 1
    param[7] = 0
    param[8] = 0
    param[9] = sp[1]
    param[10]= sp[2]
    param[11]= sp[3]
    param[12]= 1
    param[13]= sp[4]
    param[14]= 0.5

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
  } else if (method == 'heter1') {
    param <- rep(0, 18)
    param[1] = 0
    param[2] = 0
    param[3] = 0
    param[4] = 0
    param[5] = 1
    param[6] = 1
    param[7] = 0
    param[8] = 0
    param[9] = sp[1]
    param[10]= sp[2]
    param[11]= sp[3]
    param[12]= 1
    param[13]= sp[4]
    param[14]= 0.5
    param[15]= 0
    param[16]= 0
    param[17]= 0
    param[18]= sp[5]

    paramsave <- matrix(0, niter, 19)
    mod <- .Fortran("BiQRGradientH1f",
                    y = as.double(y),
                    R = as.integer(R),
                    x = as.double(x),
                    tau = as.double(tau),
                    n = as.integer(n),
                    niter = as.integer(niter),
                    param = as.double(param),
                    paramsave = as.double(paramsave)
                    )
  } else if (method == 'heter2') {
    xdim <- dim(X)[2]
    if (is.null(init)){
      require(quantreg)
      lmcoef1 <- coef(rq(y[,1] ~ X[,-1], tau = tau))
      lmcoef2 <- coef(rq(y[,2][R == 1] ~ X[R == 1,-1], tau = tau))
      param <- rep(0, 8*xdim + 3)
      param[1:xdim] <- lmcoef1
      param[(4*xdim + 1):(5*xdim)] <- lmcoef2
      param[(5*xdim + 1):(6*xdim)] = sp[1:xdim]
      param[(7*xdim + 1):(8*xdim)] = sp[(xdim + 2):(2*xdim + 1)]
      param[8*xdim + 2]= sp[xdim + 1]
      param[8*xdim + 3] = sum(R)/dim(X)[1]
    } else {
      param <- init
    }

##    print(param)
    paramsave <- matrix(0, niter, 8*xdim + 4)
    res <- matrix(0, dim(X)[1], 2)
    mod <- .Fortran("BiQRGradientH2f",
                    y = as.double(y),
                    R = as.integer(R),
                    x = as.double(X),
                    tau = as.double(tau),
                    n = as.integer(n),
                    niter = as.integer(niter),
                    param = as.double(param),
                    paramsave = as.double(paramsave),
                    xdim = as.integer(xdim),
                    converge = as.logical(TRUE),
                    res = as.double(res),
                    iter = as.integer(1),
                    tol = as.double(tol)
                    )
  }
  mod$res <- matrix(mod$res, dim(X)[1], 2)
  mod$method <- method
  class(mod) <- "BiQRGradient"
  return(mod)
}

coef.BiQRGradient <- function(mod, ...){
  q <- mod$xdim
  param <- mod$param[1:(8*q)]
  coef <- matrix(param, 8, q, byrow = T)
  coef <- coef[c(1, 5, 2, 6, 3, 4, 7, 8), ]
  p <- mod$param[8*q + 3]
  beta22 <- mod$param[8*q + 1]
  h <- mod$param[8*q + 2]
  rownames(coef) <- c('Q1Coef', 'Q2Coef', 'beta1(0)', 'beta2(0)',
                      'Sigma1(1)', 'Sigma1(0)', 'Sigma2(1)', 'Sigma2(0)')
  return(coef[c(1, 2),])
}

print.BiQRGradient <- function(mod, ...){
  cat('Coefficients: \n')
  print(coef(mod))
}

summary.BiQRGradient <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  tau <- mod$tau
  param <- mod$param
  q <- mod$xdim

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', param[8*q + 3], '\n')
  cat('Quantile: ', tau, '\n')
  cat('Quantile regression coefficients: \n')
  print(coef(mod))
}

plot.BiQRGradient <- function(mod, ...){
  a <- matrix(mod$paramsave, mod$niter, length(mod$param) + 1)
  xdim <- mod$xdim
  colnames(a) <- c(paste('gamma1', 1:xdim, sep=''), paste('beta1', 1:xdim, sep=''), paste('sigma1', 1:xdim, sep = ''), paste('sigma0', 1:xdim, sep = ''), paste('gamma2', 1:xdim, sep=''), paste('beta1', 1:xdim, sep=''), paste('sigma2', 1:xdim, sep = ''), paste('sp', 1:xdim, sep = ''), 'beta22', 'h', 'p', 'nll')
  par(mfrow = c(2, 2))
  for (i in 1:dim(a)[2]){
    plot(ts(a[1:mod$iter, i]), main = colnames(a)[i])
  }
}

Diagnose <- function(mod){
  a <- matrix(mod$paramsave, mod$niter, length(mod$param) + 1)
  xdim <- mod$xdim
  colnames(a) <- c(paste('gamma1', 1:xdim, sep=''), paste('beta1', 1:xdim, sep=''), paste('sigma1', 1:xdim, sep = ''), paste('sigma0', 1:xdim, sep = ''), paste('gamma2', 1:xdim, sep=''), paste('beta1', 1:xdim, sep=''), paste('sigma2', 1:xdim, sep = ''), paste('sp', 1:xdim, sep = ''), 'beta22', 'h', 'p', 'nll')
  par(mfrow = c(2, 2))
  for (i in 1:dim(a)[2]){
    plot(ts(a[, i]), main = colnames(a)[i])
  }
}

GoF <- function(mod){
  R <- mod$R
  res <- mod$res
  tau <- mod$tau
  qqnorm(res[, 1], main = bquote(paste('Normal Q-Q Plot for ', Y[1], ' of quantile ', tau == .(tau))))
  qqline(res[, 1])
  qqnorm(res[R == 1, 1], main = bquote(paste('Normal Q-Q Plot for ', Y[2], ' of quantile ', tau == .(tau))))
  qqline(res[R == 1, 1])
}
