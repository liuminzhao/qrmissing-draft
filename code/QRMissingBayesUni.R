#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/17/2013 17:43:47>
##' 2013/08/17 QRMissing Univariate Bayesian single normal
##' QRMissingBayesUni.f

QRMissingBayesUni <- function(y, R, X, tau = 0.5,
                              mcmc, prior
                              ){

  if (!is.loaded('qrmissingbayesuni')){
    dyn.load("~/Documents/qrmissing/code/QRMissingBayesUni.so")
  }

  ## data
  n <- length(y)
  xdim <- dim(X)[2]

  ## prior
  betapm <- prior$betapm
  betapv <- prior$betapv
  gammapm <- prior$gammapm
  gammapv <- prior$gammapv
  a <- prior$a
  b <- prior$b
  c <- prior$c
  d <- prior$d

  ## MCMC
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nburn <- mcmc$nburn
  ndisp <- mcmc$ndisp

  ## SAVE
  betasave <- gammasave <- matrix(0, nsave, xdim)
  sigmasave <- matrix(0, nsave, 2)
  psave <- rep(0, nsave)

  ## MCMC
  mod <- .Fortran("qrmissingbayesuni",
                  n = as.integer(n),
                  xdim = as.integer(xdim),
                  X = as.double(X),
                  y = as.double(y),
                  R = as.integer(R),
                  tau = as.double(tau),
                  gammapm = as.double(gammapm),
                  gammapv = as.double(gammapv),
                  betapm = as.double(betapm),
                  betapv = as.double(betapv),
                  a = as.double(a),
                  b = as.double(b),
                  c = as.double(c),
                  d = as.double(d),
                  mcmc = as.integer(mcmc),
                  nsave = as.integer(nsave),
                  gammasave = as.double(gammasave),
                  betasave = as.double(betasave),
                  sigmasave = as.double(sigmasave),
                  psave = as.double(psave)
                  )

  gammasave <- matrix(mod$gammasave, nsave, xdim)
  betasave <- matrix(mod$betasave, nsave, xdim)
  sigmasave <- matrix(mod$sigmasave, nsave, 2)
  psave <- mod$psave

  ans <- list(gammasave = gammasave,
              betasave = betasave,
              sigmasave = sigmasave,
              psave = psave,
              n = n,
              xdim = xdim,
              y = y,
              X = X,
              R = R,
              tau = tau
              )

  class(ans) <- 'QRMissingBayes'

  return(ans)

}

coef.QRMissingBayes <- function(mod, ...){
  gamma <- apply(mod$gammasave, 2, mean)
  beta <- apply(mod$betasave, 2, mean)
  sigma <- apply(mod$sigmasave, 2, mean)
  p <- mean(mod$psave)
  return(list(gamma = gamma, beta = beta, sigma = sigma, p = p))
}

summary.QRMissingBayes <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  tau <- mod$tau
  param <- mod$par
  q <- mod$xdim

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', coef(mod)$p, '\n')
  cat('Quantile: ', tau, '\n')
  cat('Quantile regression coefficients: \n')
  print(coef(mod))
}

plot.QRMissingBayes <- function(mod, ...){
  xdim <- mod$xdim
  for (i in 1:xdim){
    plot(ts(mod$gammasave[, i]), main = paste('gamma', i, sep = ''))
  }
  for (i in 1:xdim){
    plot(ts(mod$betasave[, i]), main = paste('beta', i, sep = ''))
  }
  for (i in 1:2){
    plot(ts(mod$sigmasave[, i]), main = paste('sigma', i, sep = ''))
  }
  plot(ts(mod$psave), main = 'p')
}
