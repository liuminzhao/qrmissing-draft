#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/22/2013 17:48:21>
##' 2013/08/17 QRMissing Univariate Bayesian single normal
##' QRMissingBayesUni.f
##' 2013/08/19 using pure R

ll <- function(gamma, beta, sigma, p, tau, y, X, R){
  n <- length(y)
  xdim <- dim(X)[2]
  num <- sum(R)

  d <- rep(0, n)
  d <- .Fortran("mydelta",
                x = as.double(X),
                gamma = as.double(gamma),
                beta = as.double(beta),
                sigma = as.double(sigma),
                p = as.double(p),
                tau = as.double(tau),
                n = as.integer(n),
                xdim = as.integer(xdim),
                delta = as.double(d))$delta

  lp <- X %*% beta
  mu11 <- d + lp
  mu10 <- d - lp
  ll11 <- sum(dnorm(y, mu11, sigma[1], log=T)[R==1])
  ll10 <- sum(dnorm(y, mu10, sigma[2], log=T)[R==0])
  ans <- ll11 + ll10 + num*log(p) + (n - num)*log(1 - p)

  return(ans)
}

QRMissingBayesUni <- function(y, R, X, tau = 0.5,
                              mcmc, prior
                              ){

  if (!is.loaded('mydelta')){
    dyn.load("~/Documents/qrmissing/code/QRMissingBayesUni.so")
  }

  ## data
  n <- length(y)
  xdim <- dim(X)[2]
  num <- sum(R)

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

  ## TUNE
  tunegamma <- tunebeta <- rep(0.3, xdim)
  tunesigma <- 0.3
  tunep <- 0.1
  arate <- 0.25
  attgamma <- accgamma <- attbeta <- accbeta <- rep(0, xdim)
  attsigma <- attp <- accsigma <- accp <- 0

  ## initial

  gamma <- beta <- rep(0, xdim)
  sigma <- rep(1, 2)
  p <- num/n

  arate <- 0.25

  isave <- 0
  skipcount <- 0
  dispcount <- 0
  nscan <- nburn + nskip * nsave

  start <- proc.time()[3]
  ## first

  loglikeo <- ll(gamma, beta, sigma, p, tau,  y, X, R)

  ## roll

  for (iscan in 1:nscan) {
    ## gamma
    for (i in 1:xdim){
      attgamma[i] = attgamma[i] + 1
      gammac <- gamma
      gammac[i] <- rnorm(1, gamma[i], tunegamma[i])
      logpriorc <- dnorm(gammac[i], gammapm[i], gammapv[i], log = T)
      logprioro <- dnorm(gamma[i], gammapm[i], gammapv[i], log = T)

      loglikec <- ll(gammac, beta, sigma, p, tau, y, X, R)

      ratio <- loglikec + logpriorc - loglikeo - logprioro

      if (log(runif(1)) <= ratio) {
        accgamma[i] <- accgamma[i] + 1
        loglikeo <- loglikec
        gamma <- gammac
      }
    }

    ## beta
    for (i in 1:xdim){
      attbeta[i] <- attbeta[i] + 1
      betac <- beta
      betac[i] <- rnorm(1, beta[i], tunebeta[i])
      logpriorc <- dnorm(betac[i], betapm[i], betapv[i], log = T)
      logprioro <- dnorm(beta[i], betapm[i], betapv[i], log = T)

      loglikec <- ll(gamma, betac, sigma, p, tau, y, X, R)

      ratio <- loglikec + logpriorc - loglikeo - logprioro

      if (log(runif(1)) <= ratio) {
        accbeta[i] <- accbeta[i] + 1
        loglikeo <- loglikec
        beta <- betac
      }
    }

    ## sigma
    attsigma <- attsigma + 1
    theta <- log(sigma)
    thetac <- rnorm(2, theta, tunesigma)
    logcgkc <- -sum(theta)
    logcgko <- -sum(thetac)
    sigmac <- exp(thetac)

    loglikec <- ll(gamma, beta, sigmac, p, tau, y, X, R)

    logpriorc <- sum(dgamma(sigmac, a/2, b/2, log = T))
    logprioro <- sum(dgamma(sigma, a/2, b/2, log = T))

    ratio <- loglikec + logpriorc - loglikeo - logprioro + logcgkc - logcgko

    if (log(runif(1)) <= ratio) {
      accsigma <- accsigma + 1
      loglikeo <- loglikec
      sigma <- sigmac
    }

    ## p
    attp <- attp + 1
    pc <- rnorm(1, p, tunep)
    pc <- max(min(pc, 0.99), 0.01)

    loglikec <- ll(gamma, beta, sigma, pc, tau, y, X, R)

    logpriorc <- dbeta(pc, c/2, d/2, log = T)
    logprioro <- dbeta(p, c/2, d/2, log = T)

    ratio=loglikec + logpriorc -loglikeo -logprioro

    if (log(runif(1)) <= ratio) {
      accp <- accp + 1
      loglikeo <- loglikec
      p <- pc
    }

    ## TUNE
    if (attgamma[1] >= 100 && iscan < nburn) {
      tunegamma <- tunegamma*ifelse(accgamma/attgamma > arate,
                                    2, 0.5)
      tunebeta <- tunebeta*ifelse(accbeta/attbeta > arate,
                                    2, 0.5)
      tunesigma <- tunesigma*ifelse(accsigma/attsigma > arate, 2, 0.5)
      tunep <- tunep*ifelse(accp/attp > arate, 2, 0.5)
      attgamma <- accgamma <- attbeta <- accbeta <- rep(0, xdim)
      attsigma <- accsigma <- attp <- accp <- 0
##      cat(tunegamma, tunebeta, tunesigma, tunep)
    }

    ## save

    if (iscan > nburn) {
      skipcount = skipcount + 1
      if (skipcount >= nskip) {
        isave <- isave + 1
        dispcount <- dispcount + 1
        gammasave[isave, ] <- gamma
        betasave[isave, ] <- beta
        sigmasave[isave, ] <- sigma
        psave[isave] <- p
        skipcount <- 0
        if (dispcount >= ndisp) {
          dispcount <- 0
          cat(isave, 'out of', nsave, proc.time()[3] - start, '\n')
        }
      }
    }
  }

  ans <- list(gammasave = gammasave,
              betasave = betasave,
              sigmasave = sigmasave,
              psave = psave,
              n = n,
              xdim = xdim,
              y = y,
              X = X,
              R = R,
              tau = tau,
              tune = list(gamma = tunegamma, beta = tunebeta,
                sigma = tunesigma, p = tunep)
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
