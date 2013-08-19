Delta <- function(x, gamma, beta, sigma, p){
  quan <- x %*% gamma
  lp <- x %*% beta
  d <- quan - (2*p - 1) * lp
  return(d)
}

ll <- function(gamma, beta, sigma, p, y, X, R){
  n <- length(y)
  xdim <- dim(X)[2]
  num <- sum(R)

  d <- apply(X, 1, Delta,
             gamma = gamma, beta = beta, sigma = sigma, p = p)
  lp <- X %*% as.matrix(beta)
  mu11 <- d + lp
  mu10 <- d - lp
  ll11 <- sum(dnorm(y, mu11, sigma[1], log=T)[R==1])
  ll10 <- sum(dnorm(y, mu10, sigma[2], log=T)[R==0])
  ans <- ll11 + ll10 + num*log(p) + (n - num)*log(1 - p)

  return(ans)
}

MeanMissingBayesUni <- function(y, R, X, mcmc){
  xdim <- dim(X)[2]
  n <- length(y)
  num <- sum(R)

  ## MCMC
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nburn <- mcmc$nburn
  ndisp <- mcmc$ndisp

  ## prior
  betapm <- gammapm <- rep(0, xdim)
  betapv <- gammapv <- rep(100, xdim)

  ## save
  gammasave <- betasave <- matrix(0, nsave, xdim)
  sigmasave <- matrix(0, nsave, 2)
  psave <- rep(0, nsave)

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

  loglikeo <- ll(gamma, beta, sigma, p, y, X, R)

  ## roll

  for (iscan in 1:nscan) {
    ## gamma
    for (i in 1:xdim){
      gammac <- gamma
      gammac[i] <- rnorm(1, gamma[i], 0.3)
      logpriorc <- dnorm(gammac[i], gammapm[i], gammapv[i], log = T)
      logprioro <- dnorm(gamma[i], gammapm[i], gammapv[i], log = T)

      loglikec <- ll(gammac, beta, sigma, p, y, X, R)

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

      loglikec <- ll(gamma, betac, sigma, p, y, X, R)

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

    loglikec <- ll(gamma, beta, sigmac, p, y, X, R)

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

    loglikec <- ll(gamma, beta, sigma, pc, y, X, R)

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
          cat(isave, 'out of ', nsave, '\n')
        }
      }
    }

  }

  mod <- list(gammasave = gammasave, betasave = betasave,
              sigmasave = sigmasave, psave = psave,
              n = n,
              y = y,
              xdim = xdim, X = X, R = R)

  class(mod) <- 'MeanMissingBayes'
  return(mod)
}

coef.MeanMissingBayes <- function(mod, ...){
  gamma <- apply(mod$gammasave, 2, mean)
  beta <- apply(mod$betasave, 2, mean)
  sigma <- apply(mod$sigmasave, 2, mean)
  p <- mean(mod$psave)
  return(list(gamma = gamma, beta = beta, sigma = sigma, p = p))
}

summary.MeanMissingBayes <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  q <- mod$xdim

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', coef(mod)$p, '\n')
  print(coef(mod))
}

plot.MeanMissingBayes <- function(mod, ...){
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
