#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/03/2013 22:49:11>
##' 2013/07/30 Rewrite BiMLESigma.R using pure R language
##' used uniroot.all to obtain roots
##' used optim to optimize the likelihood to get the MLE
##' used JIT/compiler package to increase speed
##' checked some small experiments and confirmed SP does not
##' affect estimates of Y1
##' 2013/08/01 try to separate D1, D2, and might run faster
##'
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Quantile Regression in the Presence of Monotone Missingness; Bivariates case
##' @param y[n, 2]
##' @param R[n]: 0-1 indicator
##' @param X[n, p]
##' @param tau: quantile requested
##' @param sp: sensitivity paramter
##' @param init: initial value for parameters
##' @param method: optimization method
##' @param tol: root finding tolerance
##' @return
##' @author Minzhao Liu
Target1 <- function(d, tau, p, quan, lp, sigma1, sigma0){
  return(tau - p*pnorm((quan - d - lp)/sigma1) - (1 - p)*pnorm(
    (quan - d + lp)/sigma0))
}

DDelta1 <- function(x, param, tau, tol, sp){
  ## translate param
  xdim <- length(x)
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

  quan <- gamma1 %*% x
  lp <- beta1 %*% x
  sigma1 <- exp(sigma11 %*% x)
  sigma0 <- exp(sigma10 %*% x)
  interval <- c(-10, 10)
  ans <- uniroot.all(Target1, interval, tau = tau, p = p ,
                     quan = quan, lp = lp,
                     sigma1 = sigma1, sigma0 = sigma0, tol = tol)[1]
  rootiter <- 0
  repeat {
    if (!is.na(ans)) {
      break
    } else {
      interval <- interval * 2
      ## cat('There is NA in D1 root \n')
      ans <- uniroot.all(Target1, interval, tau = tau, p = p ,
                         quan = quan, lp = lp,
                         sigma1 = sigma1, sigma0 = sigma0, tol = tol)[1]
    }
    rootiter <- rootiter + 1
    if (rootiter > 50) {
      cat('can not bracket root fot d1 \n')
      break
    }
  }

  ## TODO bracket low and upp until get root
  return(ans)
}

Target2 <- function(d2, x, param, tau, d1, sp){
  ## translate param
  xdim <- length(x)
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

  quan1 <- gamma1 %*% x
  lp1 <- beta1 %*% x
  sigma1 <- exp(sigma11 %*% x)
  sigma0 <- exp(sigma10 %*% x)
  mu11 <- d1 + lp1
  mu10 <- d1 - lp1
  quan2 <- gamma2 %*% x
  lp2 <- beta2 %*% x
  sigma2 <- exp(sigma21 %*% x)
  sigma2sp <- exp(sigma20 %*% x)

  if (betay == 0){
    int1 <- pnorm((quan2 - d2)/sigma2)
  } else {
    int1 <- pnorm(((quan2 - d2)/betay - mu11)/sqrt(sigma2^2/betay^2 + sigma1^2))
    if (betay < 0) {
      int1 <- 1 - int1
    }
  }

  if (beta2y == 0){
    int2 <- pnorm((quan2 - d2 - lp2)/sigma2sp)
  } else {
    int2 <- pnorm(((-d2 - lp2 + quan2)/beta2y - mu10)/sqrt(sigma2sp^2/beta2y^2 + sigma0^2))
    if (beta2y < 0){
      int2 <- 1 - int2
    }
  }
  return(tau - p*int1 - (1-p)*int2)
}

DDelta2 <- function(x, param, tau, tol, sp){
  d1 <- DDelta1(x, param, tau, tol, sp)
  interval <- c(-10, 10)
  ans <- c(d1, uniroot.all(Target2, interval,
                           x = x, param = param, tau = tau,
                           d1 = d1, sp = sp, tol = tol)[1])
  rootiter <- 0
  repeat {
    if (!is.na(ans[2])) {
      break
    } else {
      interval <- interval * 2
##           cat('THere is NA in D2 root. \n')
  ans <- c(d1, uniroot.all(Target2, interval,
                           x = x, param = param, tau = tau,
                           d1 = d1, sp = sp, tol = tol)[1])
    }
    rootiter <- rootiter + 1
    if (rootiter > 50) {
      cat('can not bracket the root for d2 \n')
      break
    }
  }
  return(ans)
}

## negative log likelihood function
nll <- function(param, y, X, R, tau, sp, tol){
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
  delta <- t(apply(X, 1, function(l) DDelta2(l, param, tau , tol, sp)))
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


QRMissingBi <- function(y, R, X, tau = 0.5, sp = NULL,
                        init = NULL, method = 'BFGS',
                        tol = 0.00001, maxit = 1000,
                        trace = 0){

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
    lmcoef1 <- coef(rq(y[,1] ~ X[,-1], tau = tau))
    lmcoef2 <- coef(rq(y[,2][R == 1] ~ X[R == 1,-1], tau = tau))
    param <- rep(0, 6*xdim + 2)
    param[1:xdim] <- lmcoef1
    param[(4*xdim + 1):(5*xdim)] <- lmcoef2
    param[6*xdim + 2] = log(num/(n-num))
  }

  nllc <- cmpfun(nll)
  ## optimize nll to get MLE
  ##  mod <- optim(param, nllc, method = method, control = list(maxit = maxit))
  mod <- optim(param, nllc, y = y, X = X, R = R, tau = tau, sp = sp, tol = tol,  method = method, control = list(maxit = maxit, trace = trace))

  mod$n <- n
  mod$xdim <- xdim
  mod$X <- X
  mod$y <- y
  mod$R <- R
  mod$tau <- tau

  class(mod) <- "QRMissingBi"

  return(mod)

}

coef.QRMissingBi <- function(mod, ...){
  q <- mod$xdim
  param <- mod$par[c(1:q, (4*q + 1):(5*q))]
  coef <- matrix(param, 2, q, byrow = T)
  rownames(coef) <- c('Q1Coef', 'Q2Coef')
  return(coef)
}

print.QRMissingBi <- function(mod, ...){
  cat('Coefficients: \n')
  print(coef(mod))
}

summary.QRMissingBi <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  tau <- mod$tau
  param <- mod$par
  q <- mod$xdim

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', exp(param[6*q + 2])/(1 + exp(param[6*q + 2])), '\n')
  cat('Quantile: ', tau, '\n')
  cat('Model converged: ', ifelse(mod$convergence, 'No', 'Yes'), '\n')
  cat('Quantile regression coefficients: \n')
  print(coef(mod))

}
