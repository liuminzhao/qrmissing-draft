#!/bin/Rscript
##' Time-stamp: <liuminzhao 07/14/2013 23:43:23>
##' 2013/07/05

ToursMNAR <- function(y, R, X, tau = 0.5, niter = 1000,  method="heter2"){
  if (!is.loaded('toursmnarh2f')) {
    if (file.exists('toursmnar.so')) {
      dyn.load('toursmnar.so')
    } else {
      if (file.exists('~/Documents/qrmissing/code/toursmnar.so')) {
        dyn.load("~/Documents/qrmissing/code/toursmnar.so")
      } else stop('no shared library found.')
    }
  }
  n <- length(R)
  xdim <- dim(X)[2]
  sp <- rep(0, 2*xdim + 1)
  sp[1] <- 3.6/100
  require(quantreg)
  lmcoef1 <- coef(rq(y[,1] ~ X[,-1], tau = tau))
  lmcoef2 <- coef(rq(y[,2][R == 1] ~ X[R == 1,-1], tau = tau))
  if (tau == 0.1) {
    lmcoef1 <- c(0.080726637783397, 0.00451080465161577, -0.442936708038249, 0.753033829630134)
    lmcoef2 <- c(-0.206850522306966, -0.00232449257256374, -0.0183695937996462, 0.310703921030932)
  }
  if (tau == 0.9) {
    lmcoef1 <- c(0.112940628486934,  0.00540403561172756, -0.0692974826323355, 0.927829207235262)
    lmcoef2 <- c(0.217754080961376, -0.0026205522008995, -0.00553671984010173, 0.325975371558977)
  }
  param <- rep(0, 8*xdim + 3)
  param[1:xdim] <- lmcoef1
  param[(4*xdim + 1):(5*xdim)] <- lmcoef2
  param[(5*xdim + 1):(6*xdim)] = sp[1:xdim]
  param[(7*xdim + 1):(8*xdim)] = sp[(xdim + 2):(2*xdim + 1)]
  param[8*xdim + 2]= sp[xdim + 1]
  param[8*xdim + 3] = 0.5
  ##    print(param)
  paramsave <- matrix(0, niter, 8*xdim + 4)
  mod <- .Fortran("toursmnarh2f",
                  y = as.double(y),
                  R = as.integer(R),
                  x = as.double(X),
                  tau = as.double(tau),
                  n = as.integer(n),
                  niter = as.integer(niter),
                  param = as.double(param),
                  paramsave = as.double(paramsave),
                  xdim = as.integer(xdim),
                  converge = as.logical(TRUE)
                  )
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
    plot(ts(a[, i]), main = colnames(a)[i])
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
