#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/19/2013 13:12:38>
##' 2013/08/17 QRMissing Univariate Bayesian single normal
##' QRMissingBayesUni.f
##' 2013/08/19 using pure R

Delta <- function(x, gamma, beta, sigma, p, tau){
  quan <- x %*% gamma
  lp <- x %*% beta
  Target <- function(d, quan, lp, sigma, p, tau){
    return(tau - p*pnorm((quan - d - lp)/sigma[1]) - (1 - p)*pnorm(
      (quan - d + lp)/sigma[2]))
  }
  interval <- c(-10, 10)
  ans <- uniroot.all(Target, interval, tau = tau, p = p ,
                     quan = quan, lp = lp,
                     sigma = sigma)[1]
  rootiter <- 0
  repeat {
    if (!is.na(ans)) {
      break
    } else {
      interval <- interval * 2
      ans <- uniroot.all(Target, interval, tau = tau, p = p ,
                         quan = quan, lp = lp,
                         sigma = sigma)[1]
    }
    rootiter <- rootiter + 1
    if (rootiter > 50) {
      print(quan)
      stop('Error')
      cat('can not bracket root fot d1 \n')
      break
    }
  }

  return(ans)
}

ll <- function(gamma, beta, sigma, p, tau, y, X, R){
  n <- length(y)
  xdim <- dim(X)[2]
  num <- sum(R)

  d <- apply(X, 1, Delta,
             gamma = gamma, beta = beta, sigma = sigma, p = p, tau = tau)
  lp <- X %*% as.matrix(beta)
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

  if (!is.loaded('qrmissingbayesuni')){
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

  ## MCMC
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
      gammac <- gamma
      gammac[i] <- rnorm(1, gamma[i], 0.3)
      logpriorc <- dnorm(gammac[i], gammapm[i], gammapv[i], log = T)
      logprioro <- dnorm(gamma[i], gammapm[i], gammapv[i], log = T)

      loglikec <- ll(gammac, beta, sigma, p, tau, y, X, R)

      ratio <- loglikec + logpriorc - loglikeo - logprioro

      if (log(runif(1)) <= ratio) {
        loglikeo <- loglikec
        gamma <- gammac
      }
    }

    ## beta
    for (i in 1:xdim){
      betac <- beta
      betac[i] <- rnorm(1, beta[i], 0.3)
      logpriorc <- dnorm(betac[i], betapm[i], betapv[i], log = T)
      logprioro <- dnorm(beta[i], betapm[i], betapv[i], log = T)

      loglikec <- ll(gamma, betac, sigma, p, tau, y, X, R)

      ratio <- loglikec + logpriorc - loglikeo - logprioro

      if (log(runif(1)) <= ratio) {
        loglikeo <- loglikec
        beta <- betac
      }
    }

    ## sigma
    theta <- log(sigma)
    thetac <- rnorm(2, theta, 0.3)
    logcgkc <- -sum(theta)
    logcgko <- -sum(thetac)
    sigmac <- exp(thetac)

    loglikec <- ll(gamma, beta, sigmac, p, tau, y, X, R)

    logpriorc <- sum(dgamma(sigmac, 1/2, 1/2, log = T))
    logprioro <- sum(dgamma(sigma, 1/2, 1/2, log = T))

    ratio <- loglikec + logpriorc - loglikeo - logprioro + logcgkc - logcgko

    if (log(runif(1)) <= ratio) {
      loglikeo <- loglikec
      sigma <- sigmac
    }

    ## p
    pc <- rnorm(1, p, 0.1)
    pc <- max(min(pc, 0.99), 0.01)

    loglikec <- ll(gamma, beta, sigma, pc, tau, y, X, R)

    logpriorc <- dbeta(pc, 1/2, 1/2, log = T)
    logprioro <- dbeta(p, 1/2, 1/2, log = T)

    ratio=loglikec + logpriorc -loglikeo -logprioro

    if (log(runif(1)) <= ratio) {
      loglikeo <- loglikec
      p <- pc
    }

    ## save

    if (iscan >= nburn) {
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

  ## mod <- .Fortran("qrmissingbayesuni",
  ##                 n = as.integer(n),
  ##                 xdim = as.integer(xdim),
  ##                 X = as.double(X),
  ##                 y = as.double(y),
  ##                 R = as.integer(R),
  ##                 tau = as.double(tau),
  ##                 gammapm = as.double(gammapm),
  ##                 gammapv = as.double(gammapv),
  ##                 betapm = as.double(betapm),
  ##                 betapv = as.double(betapv),
  ##                 a = as.double(a),
  ##                 b = as.double(b),
  ##                 c = as.double(c),
  ##                 d = as.double(d),
  ##                 mcmc = as.integer(mcmc),
  ##                 nsave = as.integer(nsave),
  ##                 gammasave = as.double(gammasave),
  ##                 betasave = as.double(betasave),
  ##                 sigmasave = as.double(sigmasave),
  ##                 psave = as.double(psave)
  ##                 )

  ## gammasave <- matrix(mod$gammasave, nsave, xdim)
  ## betasave <- matrix(mod$betasave, nsave, xdim)
  ## sigmasave <- matrix(mod$sigmasave, nsave, 2)
  ## psave <- mod$psave

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
