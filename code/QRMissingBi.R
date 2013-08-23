#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/23/2013 16:37:55>
##' 2013/07/30 Rewrite BiMLESigma.R using pure R language
##' used uniroot.all to obtain roots
##' used optim to optimize the likelihood to get the MLE
##' used JIT/compiler package to increase speed
##' checked some small experiments and confirmed SP does not
##' affect estimates of Y1
##' 2013/08/03 change rootiter from 15 to 50
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
ll2 <- function(param, y, X, R, tau, sp){
  n <- dim(y)[1]
  xdim <- dim(X)[2]
  num <- sum(R)

  gamma1 <- param[1:xdim]
  beta1 <- param[(xdim + 1):(2*xdim)]
  sigma1 <- c(exp(param[3*xdim + 3]), exp(param[3*xdim + 4]))
  gamma2 <- param[(2*xdim + 1):(3*xdim)]
  beta2sp <- sp[1:xdim] # SP for R = 0
  sigma21 <- exp(param[3*xdim + 5])
  sigma21sp <- exp(sp[xdim + 2])  # SP for R = 0
  betay <- param[3*xdim + 1] # for R = 1
  betaysp <- sp[xdim + 1] # SP for R = 0
  p <- exp(param[3*xdim + 2])/(1 + exp(param[3*xdim + 2]))

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

  return(-ans)
}


QRMissingBi <- function(y, R, X, tau = 0.5, sp = NULL,
                        init = NULL, method = 'uobyqa',
                        tol = 0.00001, control = list(maxit = 1000,
                        trace = 0), hess = FALSE){
  if (!is.loaded('mydelta2')){
    dyn.load("~/Documents/qrmissing/code/QRMissingBayesBi.so")
  }

  ## data
  n <- dim(y)[1]
  num <- sum(R)
  xdim <- dim(X)[2]

  ## initial
  if (is.null(sp)) {
    sp <- rep(0, xdim  + 2)
  }
  if (!is.null(init)){
    param <- init
  } else {
    lmcoef1 <- coef(rq(y[,1] ~ X[,-1], tau = tau))
    lmcoef2 <- coef(rq(y[,2][R == 1] ~ X[R == 1,-1], tau = tau))
    param <- rep(0, 3*xdim + 5)
    param[1:xdim] <- lmcoef1
    param[(2*xdim + 1):(3*xdim)] <- lmcoef2
    param[3*xdim + 2] = log(num/(n-num))
  }


  ## nll
  nll <- function(param){
    ll2(param, y, X, R, tau, sp)
  }

  ## optimize nll to get MLE
  optim_method <- c('BFGS', 'CG', 'L-BFGS-B', 'Nelder-Mead')

  if (method %in% optim_method) {
    mod <- optim(param, nll, method = method, control = control)
  } else {
    minqa_control <- list(iprint = control$trace)
    if (method == 'bobyqa'){
      mod <- bobyqa(param, nll, control=minqa_control)
    } else if (method == 'uobyqa') {
      mod <- uobyqa(param, nll, control=minqa_control)
    } else if (method == 'newuoa') {
      mod <- newuoa(param, nll, control=minqa_control)
    }
  }

  ## Hessian matrix and grad
  if (hess) {
    Hessian <- hessian(nll, mod$par)
    d <- grad(nll, mod$par)
    Jninv <- solve(Hessian)
    se <- matrix(0, 2, xdim)
    se[1, ] <- sqrt(diag(Jninv)[1:xdim])
    se[2, ] <- sqrt(diag(Jninv)[(2*xdim + 1):(3*xdim)])
    rownames(se) <- c('Q1', 'Q2')
  } else {
    Hessian <- NULL
    se <- NULL
  }

  mod$n <- n
  mod$xdim <- xdim
  mod$X <- X
  mod$y <- y
  mod$R <- R
  mod$tau <- tau
  mod$method <- method
  mod$Hessian <- Hessian
  mod$se <- se
##   mod$d <- d

  class(mod) <- "QRMissingBi"

  return(mod)

}

coef.QRMissingBi <- function(mod, ...){
  q <- mod$xdim
  param <- mod$par[c(1:q, (2*q + 1):(3*q))]
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
  cat('Estimated pi:', exp(param[3*q + 2])/(1 + exp(param[3*q + 2])), '\n')
  cat('Quantile: ', tau, '\n')
  cat('Optimization method: ', mod$method, '\n')
  optim_method <- c('BFGS', 'CG', 'L-BFGS-B', 'Nelder-Mead')

  if (mod$method %in% optim_method) {
    cat('Model converged: ', ifelse(mod$convergence, 'No', 'Yes'), '\n')
  } else {
    cat('Model converged: ', ifelse(mod$ierr, 'No', 'Yes'), '\n')
  }
  cat('Quantile regression coefficients: \n')
  print(coef(mod))
  cat('Standard error: \n')
  print(mod$se)

}
