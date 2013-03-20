## univariate
## missing
## general
## with covariates
## MLE with gradient descent method
## Time-stamp: <liuminzhao 03/19/2013 20:16:15>

## simulate data
## y | S = 1 ~ N(\delta + b0 + x*b1, sigma1)
## y | S = 0 ~ N(\delta - b0 - x*b1, sigma2)
## P(S = 1) = pi

## y | S = 1 ~ N(b01 + x*b11, sigma1)
## y | S = 0 ~ N(b00 + x*b10, sigma2)

library(rootSolve)
TargetEqn <- function(delta, gamma, beta, sigma, tau, p, x){
  quan <- gamma[1] + gamma[2] * x
  lp <- beta[1] + beta[2] * x
  return(tau - p * pnorm((quan - delta - lp)/sigma[1]) - (1 - p) * pnorm((quan - delta + lp)/sigma[2]))
}

SolveDelta <- function(gamma, beta, sigma, tau, p, x){
  return(uniroot.all(TargetEqn, c(-30, 30), tol = 0.0001, gamma = gamma, beta = beta, sigma = sigma, tau = tau, p = p, x = x))
}

LogLikelihood <- function(y, S, x, delta, beta, sigma){
  ll <- 0
  mu1 <- delta + beta[1] + beta[2] * x
  mu0 <- delta - beta[1] - beta[2] * x
  ll1 <- dnorm(y, mean = mu1, sd = sigma[1], log = T)
  ll0 <- dnorm(y, mean = mu0, sd = sigma[2], log = T)
  ll <- c(ll1[S == 1], ll0[S == 0])
  return(-sum(ll))
}

PartialG0 <- function(gamma, beta){
  epsilon <- 0.01
  gamma1 <- c(gamma[1] + epsilon, gamma[2])
  gamma2 <- c(gamma[1] - epsilon, gamma[2])
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma1, beta, sigma, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma2, beta, sigma, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta, sigma)
  ll2 <- LogLikelihood(y, S, x, delta2, beta, sigma)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialG1 <- function(gamma, beta){
  epsilon <- 0.01
  gamma1 <- c(gamma[1], gamma[2] + epsilon)
  gamma2 <- c(gamma[1], gamma[2] - epsilon)
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma1, beta, sigma, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma2, beta, sigma, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta, sigma)
  ll2 <- LogLikelihood(y, S, x, delta2, beta, sigma)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialB0 <- function(gamma, beta){
  epsilon <- 0.01
  beta1 <- c(beta[1] + epsilon, beta[2])
  beta2 <- c(beta[1] - epsilon, beta[2])
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta1, sigma, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta2, sigma, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta1, sigma)
  ll2 <- LogLikelihood(y, S, x, delta2, beta2, sigma)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialB1 <- function(gamma, beta){
  epsilon <- 0.01
  beta1 <- c(beta[1], beta[2] + epsilon)
  beta2 <- c(beta[1], beta[2] - epsilon)
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta1, sigma, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta2, sigma, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta1, sigma)
  ll2 <- LogLikelihood(y, S, x, delta2, beta2, sigma)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialS0 <- function(gamma, beta, sigma){
  epsilon <- 0.001
  sigma1 <- c(sigma[1] + epsilon, sigma[2])
  sigma2 <- c(sigma[1] - epsilon, sigma[2])
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, sigma1, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, sigma2, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta, sigma1)
  ll2 <- LogLikelihood(y, S, x, delta2, beta, sigma2)
  return((ll1 - ll2) / 2/ epsilon)
}

PartialS1 <- function(gamma, beta, sigma){
  epsilon <- 0.001
  sigma1 <- c(sigma[1], sigma[2] + epsilon)
  sigma2 <- c(sigma[1], sigma[2] - epsilon)
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, sigma1, tau, phi, x)))
  delta2 <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, sigma2, tau, phi, x)))
  ll1 <- LogLikelihood(y, S, x, delta1, beta, sigma1)
  ll2 <- LogLikelihood(y, S, x, delta2, beta, sigma2)
  return((ll1 - ll2) / 2/ epsilon)
}

QRGradientUni <- function(y, S, x, tau){
  phi <- 0.5
  beta <- c(0, 0) # beta^(1)
  gamma <- c(0, 0)
  sigma <- c(1, 1)
  delta <- rep(0, n)
  ll0 <- LogLikelihood(y, S, x, delta, beta, sigma)
  dif <- 1
  alpha <- 0.003
  iter <- 0
  llsave <- ll0
  gamma0save <- gamma1save <- beta0save <- beta1save <- sigma1save <- sigma0save <- 0
###############
  ## BEGIN GRADIENT DESCENT 
###############
  while (dif > 10^-5 & iter < 1000) {
    pg0 <- PartialG0(gamma, beta)
    pg1 <- PartialG1(gamma, beta)
    pb0 <- PartialB0(gamma, beta)
    pb1 <- PartialB1(gamma, beta)
    ps0 <- PartialS0(gamma, beta, sigma)
    ps1 <- PartialS1(gamma, beta, sigma)
    gamma[1] <- gamma[1] - alpha * pg0
    gamma[2] <- gamma[2] - alpha * pg1
    beta[1] <- beta[1] - alpha * pb0
    beta[2] <- beta[2] - alpha * pb1
    sigma[1] <- max(sigma[1] - alpha * ps0, 0.01)
    sigma[2] <- max(sigma[2] - alpha * ps1, 0.01)
    delta <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, sigma, tau, phi, x)))
    ll <- LogLikelihood(y, S, x, delta, beta, sigma)
    dif <- abs(ll - ll0)
    ll0 <- ll
    llsave <- append(llsave, ll)
    gamma0save <- append(gamma0save, gamma[1])
    gamma1save <- append(gamma1save, gamma[2])
    beta0save <- append(beta0save, beta[1])
    beta1save <- append(beta1save, beta[2])
    sigma0save <- append(sigma0save, sigma[1])
    sigma1save <- append(sigma1save, sigma[2])
    iter <- iter + 1
  }
  return(list(gamma = gamma))
}
  
