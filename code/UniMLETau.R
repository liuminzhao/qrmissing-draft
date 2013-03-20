## Time-stamp: <liuminzhao 03/20/2013 13:48:40>
## Univariate
## Using tau instead of sigma
## MLE with gradient descent
library(rootSolve)
TargetEqn <- function(delta, gamma, beta, prec, tau, p, x){
  quan <- gamma[1] + gamma[2] * x
  lp <- beta[1] + beta[2] * x
  return(tau - p * pnorm((quan - delta - lp) * sqrt(prec[1])) - (1 - p) * pnorm((quan - delta + lp) * sqrt(prec[2])))
}

SolveDelta <- function(gamma, beta, prec, tau, p, x){
  return(uniroot.all(TargetEqn, c(-30, 30), tol = 0.0001, gamma = gamma, beta = beta, prec = prec, tau = tau, p = p, x = x)[1])
}

LogLikelihood <- function(y, S, x, delta, beta, prec){
  ll <- 0
  mu1 <- delta + beta[1] + beta[2] * x
  mu0 <- delta - beta[1] - beta[2] * x
  ll1 <- dnorm(y, mean = mu1, sd = 1/sqrt(prec[1]), log = T)
  ll0 <- dnorm(y, mean = mu0, sd = 1/sqrt(prec[2]), log = T)
  ll <- c(ll1[S == 1], ll0[S == 0])
  return(-sum(ll))
}

PartialG0 <- function(gamma, beta, prec, tau, phi, y, S, x){
  epsilon <- 0.01
  gamma1 <- c(gamma[1] + epsilon, gamma[2])
  gamma2 <- c(gamma[1] - epsilon, gamma[2])
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma1, beta, prec, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma2, beta, prec, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta, prec)
  ll2 <- LogLikelihood(y, S, x, delta2, beta, prec)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialG1 <- function(gamma, beta, prec, tau, phi, y, S, x){
  epsilon <- 0.01
  gamma1 <- c(gamma[1], gamma[2] + epsilon)
  gamma2 <- c(gamma[1], gamma[2] - epsilon)
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma1, beta, prec, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma2, beta, prec, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta, prec)
  ll2 <- LogLikelihood(y, S, x, delta2, beta, prec)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialB0 <- function(gamma, beta, prec, tau, phi, y, S, x){
  epsilon <- 0.01
  beta1 <- c(beta[1] + epsilon, beta[2])
  beta2 <- c(beta[1] - epsilon, beta[2])
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta1, prec, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta2, prec, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta1, prec)
  ll2 <- LogLikelihood(y, S, x, delta2, beta2, prec)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialB1 <- function(gamma, beta, prec, tau, phi, y, S, x){
  epsilon <- 0.01
  beta1 <- c(beta[1], beta[2] + epsilon)
  beta2 <- c(beta[1], beta[2] - epsilon)
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta1, prec, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta2, prec, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta1, prec)
  ll2 <- LogLikelihood(y, S, x, delta2, beta2, prec)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialP0 <- function(gamma, beta, prec, tau, phi, y, S, x){
  epsilon <- 0.001
  prec1 <- c(prec[1] + epsilon, prec[2])
  prec2 <- c(prec[1] - epsilon, prec[2])
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, prec1, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, prec2, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta, prec1)
  ll2 <- LogLikelihood(y, S, x, delta2, beta, prec2)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialP1 <- function(gamma, beta, prec, tau, phi, y, S, x){
  epsilon <- 0.001
  prec1 <- c(prec[1], prec[2] + epsilon)
  prec2 <- c(prec[1], prec[2] - epsilon)
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, prec1, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, prec2, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta, prec1)
  ll2 <- LogLikelihood(y, S, x, delta2, beta, prec2)
  return((ll1 - ll2) / 2/ epsilon)
}

QRGradient <- function(y, S, x, tau, alpha = 0.003){

  phi <- 0.5
  n <- length(y)
  beta <- c(0, 0) # beta^(1)
  gamma <- c(0, 0)
  prec <- c(1, 1)
  delta <- rep(0, n)
  ll0 <- LogLikelihood(y, S, x, delta, beta, prec)
  dif <- 1
  alpha <- 0.003
  iter <- 0
  llsave <- ll0
  gamma0save <- gamma1save <- beta0save <- beta1save <- prec1save <- prec0save <- 0
###############
  ## BEGIN GRADIENT DESCENT 
###############
  while (dif > 10^-5 & iter < 1000) {
    pg0 <- PartialG0(gamma, beta, prec, tau, phi, y, S, x)
    pg1 <- PartialG1(gamma, beta, prec, tau, phi, y, S, x)
    pb0 <- PartialB0(gamma, beta, prec, tau, phi, y, S, x)
    pb1 <- PartialB1(gamma, beta, prec, tau, phi, y, S, x)
    pp0 <- PartialP0(gamma, beta, prec, tau, phi, y, S, x)
    pp1 <- PartialP1(gamma, beta, prec, tau, phi, y, S, x)
    gamma[1] <- gamma[1] - alpha * pg0
    gamma[2] <- gamma[2] - alpha * pg1
    beta[1] <- beta[1] - alpha * pb0
    beta[2] <- beta[2] - alpha * pb1
    prec[1] <- max(prec[1] - alpha * pp0, 0.01)
    prec[2] <- max(prec[2] - alpha * pp1, 0.01)
    delta <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, prec, tau, phi, x)))
    ll <- LogLikelihood(y, S, x, delta, beta, prec)
    dif <- abs(ll - ll0)
    ll0 <- ll
    llsave <- append(llsave, ll)
    gamma0save <- append(gamma0save, gamma[1])
    gamma1save <- append(gamma1save, gamma[2])
    beta0save <- append(beta0save, beta[1])
    beta1save <- append(beta1save, beta[2])
    prec0save <- append(prec0save, prec[1])
    prec1save <- append(prec1save, prec[2])
    iter <- iter + 1
  }

  return(list(gamma = gamma, beta = beta, prec = prec, llsave = llsave, iter = iter))
}
