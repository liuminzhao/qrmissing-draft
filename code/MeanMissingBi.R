#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/04/2013 15:23:37>
##' 2013/08/04 Mean regression for Bivariate using our approach
##' Test on SP effect on Y1 marignal estimates for mean
Delta <- function(x, param, sp){
  ## translate param
  xdim <- length(x)
  gamma1 <- param[1:xdim]
  beta1 <- param[(xdim + 1):(2*xdim)]
  gamma2 <- param[(4*xdim + 1):(5*xdim)]
  beta2 <- sp[1:xdim] # SP for R = 0
  betay <- param[6*xdim + 1] # for R = 1
  beta2y <- betay + sp[xdim + 1] # SP for R = 0
  p <- exp(param[6*xdim + 2])/(1 + exp(param[6*xdim + 2]))

  quan <- gamma1 %*% x
  lp <- beta1 %*% x

  ## d1
  d1 <- quan - (2*p - 1)*lp

  ## d2
  quan2 <- gamma2 %*% x
  lp2 <- beta2 %*% x

  d2 <- quan2 - (p*betay + (1 - p)*beta2y)*d1 -
    (p*betay - (1 - p)*beta2y)*lp - (1 - p)*lp2

  return(c(d1, d2))
}

## negative log likelihood function
nll <- function(param, y, X, R, sp){
  ## translate param
  n <- dim(y)[1]
  num <- sum(R)
  xdim <- dim(X)[2]
  xdim <- dim(X)[2]
  gamma1 <- param[1:xdim]
  beta1 <- param[(xdim + 1):(2*xdim)]
  sigma11 <- param[(2*xdim + 1):(3*xdim)]
  sigma10 <- param[(3*xdim + 1):(4*xdim)]
  gamma2 <- param[(4*xdim + 1):(5*xdim)]
  beta2 <- sp[1:xdim] # SP for R = 0
  sigma21 <- param[(5*xdim + 1):(6*xdim)]
  sigma20 <- sp[(xdim + 2):(2*xdim+1)] + sigma21 # SP for R = 0
  betay <- param[6*xdim + 1] # for R = 1
  beta2y <- betay + sp[xdim + 1] # SP for R = 0
  p <- exp(param[6*xdim + 2])/(1 + exp(param[6*xdim + 2]))

  ## get delta
  delta <- t(apply(X, 1, function(l) Delta(l, param, sp)))
  d1 <- delta[,1]
  d2 <- delta[,2]

  ## Y1
  lp1 <- X %*% as.matrix(beta1)
  mu11 <- d1 + lp1
  mu10 <- d1 - lp1
  sigma1 <- exp(X %*% as.matrix(sigma11))
  sigma0 <- exp(X %*% as.matrix(sigma10))
  ## Y2
  mu21 <- d2 + betay * y[,1]
  sigma2 <- exp(X %*% as.matrix(sigma21))

  ## ll
  ll11 <- sum(dnorm(y[,1], mu11, sigma1, log=T)[R==1])
  ll10 <- sum(dnorm(y[,1], mu10, sigma0, log=T)[R==0])
  ll21 <- sum(dnorm(y[,2], mu21, sigma2, log=T)[R==1])
  ans <- ll11+ll10+ll21+ num*log(p) + (n - num)*log(1 - p)

  return(-ans)
}


MeanMissingBi <- function(y, R, X, sp = NULL,
                          init = NULL, method = 'BFGS',
                          maxit = 1000,
                          trace = 0, lower = NULL, upper = NULL){
  ## data
  n <- dim(y)[1]
  num <- sum(R)
  xdim <- dim(X)[2]

  ## initial
  if (is.null(sp)) {
    sp <- rep(0, xdim * 2 + 1)
  }
  if (!is.null(init)){
    param <- init
  } else {
    lmcoef1 <- coef(lm(y[,1]~X[, -1]))
    lmcoef2 <- coef(lm(y[R==1,1]~X[R==1, -1]))
    param <- rep(0, 6*xdim + 2)
    param[1:xdim] <- lmcoef1
    param[(4*xdim + 1):(5*xdim)] <- lmcoef2
    param[6*xdim + 2] = log(num/(n-num))
  }

  nllc <- cmpfun(nll)
  if (is.null(lower)) {
    mod <- optim(param, nllc, y = y, X = X, R = R, sp = sp,  method = method, control = list(maxit = maxit, trace = trace))
  } else {
    mod <- optim(param, nllc, y = y, X = X, R = R, sp = sp, method = method, upper = upper, lower = lower,  control = list(maxit = maxit, trace = trace))
  }

  mod$n <- n
  mod$xdim <- xdim
  mod$X <- X
  mod$y <- y
  mod$R <- R

  class(mod) <- "MeanMissingBi"

  return(mod)

}


coef.MeanMissingBi <- function(mod, ...){
  q <- mod$xdim
  param <- mod$par[c(1:q, (4*q + 1):(5*q))]
  coef <- matrix(param, 2, q, byrow = T)
  rownames(coef) <- c('Q1Coef', 'Q2Coef')
  return(coef)
}

print.MeanMissingBi <- function(mod, ...){
  cat('Coefficients: \n')
  print(coef(mod))
}

summary.MeanMissingBi <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  param <- mod$par
  q <- mod$xdim

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', exp(param[6*q + 2])/(1 + exp(param[6*q + 2])), '\n')
  cat('Model converged: ', ifelse(mod$convergence, 'No', 'Yes'), '\n')
  cat('Quantile regression coefficients: \n')
  print(coef(mod))

}
