#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/22/2013 22:13:31>
##' 2013/08/22 QRMissing Bivariate Bayesian single normal
##' QRMissingBayesBi.f

ll2 <- function(gamma1, beta1, sigma1,
                gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
                p, tau, y, X, R){
  n <- dim(y)[1]
  xdim <- dim(X)[2]
  num <- sum(R)

  d <- matrix(0, n, 2)
  d <- .Fortran("mydelta2",
                x = as.double(X),
                gamma1 = as.double(gamma1),
                beta1 = as.double(beta1),
                sigma1 = as.double(sigma1),
                gamma2 = as.double(gamma2),
                beta2sp = as.double(beta2sp),
                sigma21 = as.double(sigma21),
                sigma21sp = as.double(sigma21sp),
                betay = as.double(betay),
                betaysp = as.double(betaysp),
                p = as.double(p),
                tau = as.double(tau),
                n = as.integer(n),
                xdim = as.integer(xdim),
                delta = as.double(d))$delta

  d <- matrix(d, n, 2)

  lp1 <- X %*% beta1
  mu11 <- d[, 1] + lp1
  mu10 <- d[, 1] - lp1
  mu21 <- d[, 2] + betay * y[, 1]
  ll11 <- sum(dnorm(y[, 1], mu11, sigma1[1], log=T)[R==1])
  ll10 <- sum(dnorm(y[, 1], mu10, sigma1[2], log=T)[R==0])
  ll21 <- sum(dnorm(y[, 2], mu21, sigma21, log = T)[R==1])
  ans <- ll11 + ll10 + ll21 + num*log(p) + (n - num)*log(1 - p)

  return(ans)
}

QRMissingBayesBi <- function(y, R, X, tau = 0.5,
                             mcmc, prior
                             ){
  require(quantreg)
  if (!is.loaded('mydelta2')){
    dyn.load("~/Documents/qrmissing/code/QRMissingBayesBi.so")
  }

  ## data
  n <- dim(y)[1]
  xdim <- dim(X)[2]
  num <- sum(R)

  ## prior
  gammapm <- prior$gammapm
  gammapv <- prior$gammapv
  betapm <- prior$betapm
  betapv <- prior$betapv
  a <- prior$a
  b <- prior$b
  c <- prior$c
  d <- prior$d

  ## SP prior
  beta2pm <- prior$beta2pm
  beta2pv <- prior$beta2pv
  sigma21sppm <- prior$sigma21sppm
  sigma21sppv <- prior$sigma21sppv
  betaysppm <- prior$betaysppm
  betaysppv <- prior$betaysppv

  ## MCMC
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nburn <- mcmc$nburn
  ndisp <- mcmc$ndisp

  ## SAVE
  beta1save <- gamma1save <- gamma2save <- matrix(0, nsave, xdim)
  sigma1save <- matrix(0, nsave, 2)
  sigma21save <- rep(0, nsave)
  betaysave <- rep(0, nsave)
  psave <- rep(0, nsave)

  ## TUNE
  tunegamma1 <- tunebeta1 <- tunegamma2 <- rep(0.3, xdim)
  tunesigma1 <- tunesigma21 <- 0.3
  tunebetay <- 0.3
  tunep <- 0.1
  arate <- 0.25
  attgamma1 <- accgamma1 <- attgamma2 <- accgamma2 <- attbeta1 <- accbeta1 <- rep(0, xdim)
  attsigma1 <- attp <- accsigma1 <- accp <- attsigma21 <- accsigma21 <- 0
  attbetay <- accbetay <- 0

  ## initial

  beta1 <- beta2 <- rep(0, xdim)
  gamma1 <- coef(rq(y[, 1] ~ X[, -1], tau = tau))
  gamma2 <- coef(rq(y[R==1,2] ~ X[R==1, -1], tau = tau))
  sigma1 <- rep(1, 2)
  sigma21 <- 1
  sigma21sp <- 0
  beta2sp <- rep(0, xdim)
  betay <- betaysp <- 0
  p <- num/n

  arate <- 0.25

  isave <- 0
  skipcount <- 0
  dispcount <- 0
  nscan <- nburn + nskip * nsave

  start <- proc.time()[3]
  ## first

  loglikeo <- ll2(gamma1, beta1, sigma1,
                  gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
                  p, tau, y, X, R)

  ## roll

  for (iscan in 1:nscan) {
    ## gamma1
    for (i in 1:xdim){
      attgamma1[i] = attgamma1[i] + 1
      gamma1c <- gamma1
      gamma1c[i] <- rnorm(1, gamma1[i], tunegamma1[i])
      logpriorc <- dnorm(gamma1c[i], gammapm[i], gammapv[i], log = T)
      logprioro <- dnorm(gamma1[i], gammapm[i], gammapv[i], log = T)

      loglikec <- ll2(gamma1c, beta1, sigma1,
                      gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
                      p, tau, y, X, R)

      ratio <- loglikec + logpriorc - loglikeo - logprioro

      ## cat(loglikeo, loglikec, gamma1, gamma1c, exp(ratio), '\n')
      ## readline("Press <return to continue")

      if (log(runif(1)) <= ratio) {
        accgamma1[i] <- accgamma1[i] + 1
        loglikeo <- loglikec
        gamma1 <- gamma1c
      }
    }

    ## beta1
    for (i in 1:xdim){
      attbeta1[i] <- attbeta1[i] + 1
      beta1c <- beta1
      beta1c[i] <- rnorm(1, beta1[i], tunebeta1[i])
      logpriorc <- dnorm(beta1c[i], betapm[i], betapv[i], log = T)
      logprioro <- dnorm(beta1[i], betapm[i], betapv[i], log = T)

      loglikec <- ll2(gamma1, beta1c, sigma1,
                      gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
                      p, tau, y, X, R)

      ratio <- loglikec + logpriorc - loglikeo - logprioro

      if (log(runif(1)) <= ratio) {
        accbeta1[i] <- accbeta1[i] + 1
        loglikeo <- loglikec
        beta1 <- beta1c
      }
    }

    ## sigma1
    attsigma1 <- attsigma1 + 1
    theta <- log(sigma1)
    thetac <- rnorm(2, theta, tunesigma1)
    logcgkc <- -sum(theta)
    logcgko <- -sum(thetac)
    sigma1c <- exp(thetac)

    loglikec <- ll2(gamma1, beta1, sigma1c,
                    gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
                    p, tau, y, X, R)


    logpriorc <- sum(dgamma(sigma1c, a/2, b/2, log = T))
    logprioro <- sum(dgamma(sigma1, a/2, b/2, log = T))

    ratio <- loglikec + logpriorc - loglikeo - logprioro + logcgkc - logcgko

    if (log(runif(1)) <= ratio) {
      accsigma1 <- accsigma1 + 1
      loglikeo <- loglikec
      sigma1 <- sigma1c
    }

    ## gamma2
    for (i in 1:xdim){
      attgamma2[i] = attgamma2[i] + 1
      gamma2c <- gamma2
      gamma2c[i] <- rnorm(1, gamma2[i], tunegamma2[i])
      logpriorc <- dnorm(gamma2c[i], gammapm[i], gammapv[i], log = T)
      logprioro <- dnorm(gamma2[i], gammapm[i], gammapv[i], log = T)

      loglikec <- ll2(gamma1, beta1, sigma1,
                      gamma2c, beta2sp, sigma21, sigma21sp, betay, betaysp,
                      p, tau, y, X, R)

      ratio <- loglikec + logpriorc - loglikeo - logprioro

      if (log(runif(1)) <= ratio) {
        accgamma2[i] <- accgamma2[i] + 1
        loglikeo <- loglikec
        gamma2 <- gamma2c
      }
    }

    ## sigma21
    attsigma21 <- attsigma21 + 1
    theta <- log(sigma21)
    thetac <- rnorm(1, theta, tunesigma21)
    logcgkc <- -theta
    logcgko <- -thetac
    sigma21c <- exp(thetac)

    loglikec <- ll2(gamma1, beta1, sigma1,
                    gamma2, beta2sp, sigma21c, sigma21sp, betay, betaysp,
                    p, tau, y, X, R)


    logpriorc <- dgamma(sigma21c, a/2, b/2, log = T)
    logprioro <- dgamma(sigma21, a/2, b/2, log = T)

    ratio <- loglikec + logpriorc - loglikeo - logprioro + logcgkc - logcgko

    if (log(runif(1)) <= ratio) {
      accsigma21 <- accsigma21 + 1
      loglikeo <- loglikec
      sigma21 <- sigma21c
    }


    ## betay
    attbetay <- attbetay + 1
    betayc <- betay
    betayc <- rnorm(1, betay, tunebetay)
    logpriorc <- dnorm(betayc, betapm[1], betapv[1], log = T)
    logprioro <- dnorm(betay, betapm[1], betapv[1], log = T)

    loglikec <- ll2(gamma1, beta1, sigma1,
                    gamma2, beta2sp, sigma21, sigma21sp, betayc, betaysp,
                    p, tau, y, X, R)

    ratio <- loglikec + logpriorc - loglikeo - logprioro

    if (log(runif(1)) <= ratio) {
      accbetay <- accbetay + 1
      loglikeo <- loglikec
      betay <- betayc
    }

    ## p
    attp <- attp + 1
    pc <- rnorm(1, p, tunep)
    pc <- max(min(pc, 0.99), 0.01)

    loglikec <- ll2(gamma1, beta1, sigma1,
                    gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
                    pc, tau, y, X, R)

    logpriorc <- dbeta(pc, c/2, d/2, log = T)
    logprioro <- dbeta(p, c/2, d/2, log = T)

    ratio=loglikec + logpriorc -loglikeo -logprioro

    if (log(runif(1)) <= ratio) {
      accp <- accp + 1
      loglikeo <- loglikec
      p <- pc
    }

    ## SP
    beta2sp <- rnorm(xdim, beta2pm, beta2pv)
    sigma21sp <- rnorm(1, sigma21sppm, sigma21sppv)
    betaysp <- rnorm(1, betaysppm, betaysppv)

    ## TUNE
    if (attgamma1[1] >= 100 & iscan < nburn) {
      tunegamma1 <- pmin(pmax(tunegamma1*ifelse(accgamma1/attgamma1 > arate,
                                                1.5, 0.7), 0.01), 10)
      tunebeta1 <- pmin(pmax(tunebeta1*ifelse(accbeta1/attbeta1 > arate,
                                    1.5, 0.7), 0.01), 10)
      tunesigma1 <- pmin(pmax(tunesigma1*ifelse(accsigma1/attsigma1 > arate, 1.5, 0.7), 0.01), 10)
      tunegamma2 <- pmin(pmax(tunegamma2*ifelse(accgamma2/attgamma2 > arate,
                                                1.5, 0.7), 0.01), 10)
      tunesigma21 <- pmin(pmax(tunesigma21*ifelse(accsigma21/attsigma21 > arate, 1.5, 0.7), 0.01), 10)
      tunebetay <- pmin(pmax(tunebetay*ifelse(accbetay/attbetay > arate,
                                              1.5, 0.7), 0.01), 10)
      tunep <- pmin(pmax(tunep*ifelse(accp/attp > arate, 1.5, 0.7), 0.01), 10)
      ## cat(accgamma1, accbeta1, accsigma1, accgamma2, accsigma21, accbetay, accp, '\n')
      ## cat(attgamma1, attbeta1, attsigma1, attgamma2, attsigma21, attbetay, attp, '\n')
      ## cat(tunegamma1, tunebeta1, tunesigma1, tunegamma2, tunesigma21, tunebetay, tunep, '\n')
      attgamma1 <- accgamma1 <- attbeta1 <- accbeta1 <- rep(0, xdim)
      attgamma2 <- accgamma2 <- rep(0, xdim)
      attsigma1 <- accsigma1 <- attsigma21 <- accsigma21 <- attp <- accp <- 0
      attsigma21 <- accsigma21 <- 0
      attbetay <- accbetay <- 0
    }

    ## save

    if (iscan > nburn) {
      skipcount = skipcount + 1
      if (skipcount >= nskip) {
        isave <- isave + 1
        dispcount <- dispcount + 1
        gamma1save[isave, ] <- gamma1
        beta1save[isave, ] <- beta1
        sigma1save[isave, ] <- sigma1
        gamma2save[isave, ] <- gamma2
        betaysave[isave] <- betay
        sigma21save[isave] <- sigma21
        psave[isave] <- p
        skipcount <- 0
        if (dispcount >= ndisp) {
          dispcount <- 0
          cat(isave, 'out of', nsave, proc.time()[3] - start, '\n')
        }
      }
    }
  }

  ans <- list(gamma1save = gamma1save,
              beta1save = beta1save,
              sigma1save = sigma1save,
              gamma2save = gamma2save,
              betaysave = betaysave,
              sigma21save = sigma21save,
              psave = psave,
              n = n,
              xdim = xdim,
              y = y,
              X = X,
              R = R,
              tau = tau,
              tune = list(gamma1 = tunegamma1, beta1 = tunebeta1,
                sigma1 = tunesigma1, p = tunep, gamma2 =tunegamma2,
                betay = tunebetay, sigma21 = tunesigma21)
              )

  class(ans) <- 'QRMissingBayesBi'

  return(ans)

}

coef.QRMissingBayesBi <- function(mod, ...){
  gamma1 <- apply(mod$gamma1save, 2, mean)
  beta1 <- apply(mod$beta1save, 2, mean)
  sigma1 <- apply(mod$sigma1save, 2, mean)
  gamma2 <- apply(mod$gamma2save, 2, mean)
  betay <- mean(mod$betaysave)
  sigma21 <- mean(mod$sigma21save)
  p <- mean(mod$psave)
  return(list(gamma1 = gamma1, beta1 = beta1, sigma1 = sigma1,
              gamma2 = gamma2, betay = betay, sigma21 = sigma21, p = p))
}

summary.QRMissingBayesBi <- function(mod, ...){
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

plot.QRMissingBayesBi <- function(mod, ...){
  xdim <- mod$xdim
  for (i in 1:xdim){
    plot(ts(mod$gamma1save[, i]), main = paste('gamma1', i, sep = ''))
  }
  for (i in 1:xdim){
    plot(ts(mod$gamma2save[, i]), main = paste('gamma2', i, sep = ''))
  }
  for (i in 1:xdim){
    plot(ts(mod$beta1save[, i]), main = paste('beta1', i, sep = ''))
  }
  for (i in 1:2){
    plot(ts(mod$sigma1save[, i]), main = paste('sigma1', i, sep = ''))
  }
  plot(ts(mod$sigma21save), main = 'Sigma21')
  plot(ts(mod$psave), main = 'p')
}
