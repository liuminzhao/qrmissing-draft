## Time-stamp: <liuminzhao 03/25/2013 20:19:58>
## MLE with Bivariate Case
## ONE COVARIATE
## FIX SENSITIVITY PARAMETER (MAR OR MNAR)

## Param
## gamma1, gamma2, beta1, beta2 (3) , prec1, prec2, lambda, phi


## FUNCTION

library(rootSolve)

## Delta1
## beta1 is beta^(1)
TargetEqn1 <- function(delta1, gamma1, beta1, prec1, tau, p, x){
  quan <- gamma1[1] + gamma1[2] * x
  lp <- beta1[1] + beta1[2] * x
  return(tau - p * pnorm((quan - delta1 - lp) * sqrt(prec1[1])) - (1 - p) * pnorm((quan - delta1 + lp) * sqrt(prec1[2])))
}

SolveDelta1 <- function(gamma1, beta1, prec1, tau, p, x){
  return(uniroot.all(TargetEqn1, c(-30, 30), tol = 0.0001, gamma1 = gamma1, beta1 = beta1, prec1 = prec1, tau = tau, p = p, x = x)[1])
}

## Delta2
## beta2 is beta_2^(0) the sensitivity parameter, prec2 = 1/ sigma2|1 ^2
TargetEqn2 <- function(delta2, gamma2, beta1, beta2, prec1, prec2, lambda, tau, p, x, delta1){
  if (beta2[3] != 0) {
    p1 <- 1 - pnorm(((delta2 - gamma2[1] - x*gamma2[2] - beta2[1] - x*beta2[2])/beta2[3] - (delta1 + beta1[1] + x*beta1[2])) * sqrt(prec1[1]) / sqrt(prec1[1]/prec2/beta2[3]^2 + 1))
    p2 <- pnorm(((-delta2 + gamma2[1] + x*gamma2[2] - beta2[1] - x*beta2[2])/beta2[3] - (delta1 - beta1[1] - x*beta1[2])) * sqrt(prec1[2]) / sqrt(prec1[2]/prec2*lambda^2/beta2[3]^2 + 1))
    if (beta2[3] > 0) {
      return(tau - p * p1 - (1 - p) * p2)
    } else {
      return(tau - p * (1 - p1) - (1- p) * (1 - p2))
    }
  } else {
    p1 <- pnorm((gamma2[1] + x*gamma2[2] - delta2 + beta2[1] + x*beta2[2]) * sqrt(prec2[1]))
    p2 <- pnorm((gamma2[1] + x*gamma2[2] - delta2 - beta2[1] - x*beta2[2]) * sqrt(prec2[1]) / lambda)
    return(tau - p * p1 - (1 - p) * p2)
  }
}

SolveDelta2 <- function(gamma2, beta1, beta2, prec1, prec2, lambda, tau, p, x, delta1){
  return(uniroot.all(TargetEqn2, c(-30, 30), tol = 0.0001, gamma2 = gamma2, beta1 = beta1, beta2 = beta2, prec1 = prec1, prec2 = prec2, lambda = lambda, tau = tau, p = p, x = x, delta1 = delta1)[1])
}

## dim(y) = c(n, 2)
NegLogLikelihood <- function(y, R, x, delta1, beta1, prec1, delta2, beta2, prec2, lambda){
  ll <- 0
  mu11 <- delta1 + beta1[1] + x*beta1[2]
  mu01 <- delta1 - beta1[1] - x*beta1[2]
  mu12 <- delta2 - beta2[1] - x*beta2[2] - y[, 1]*beta2[3]
  ll1 <- dnorm(y[, 1], mean = mu11, sd = 1/sqrt(prec1[1]), log = T)
  ll0 <- dnorm(y[, 1], mean = mu01, sd = 1/sqrt(prec1[2]), log = T)
  ll2 <- dnorm(y[, 2], mean = mu12, sd = 1/sqrt(prec2), log = T)
  ll <- c(ll1[R == 1], ll0[R == 0], ll2[R == 1])
  return(-sum(ll))
}

PartialG01 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  gamma1p <- c(gamma1[1] + epsilon, gamma1[2])
  gamma1m <- c(gamma1[1] - epsilon, gamma1[2])
  delta1p <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1p, beta1, prec1, tau, p, x)))
  delta1m <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1m, beta1, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1p)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1m)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1p, beta1, prec1, delta2p, beta2, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1m, beta1, prec1, delta2m, beta2, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialG11 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  gamma1p <- c(gamma1[1], gamma1[2] + epsilon)
  gamma1m <- c(gamma1[1], gamma1[2] - epsilon)
  delta1p <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1p, beta1, prec1, tau, p, x)))
  delta1m <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1m, beta1, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1p)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1m)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1p, beta1, prec1, delta2p, beta2, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1m, beta1, prec1, delta2m, beta2, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialG02 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  gamma2p <- c(gamma2[1] + epsilon, gamma2[2])
  gamma2m <- c(gamma2[1] - epsilon, gamma2[2])
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2p, beta1, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2m, beta1, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2p
                          , beta2, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2m
                          , beta2, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}


PartialG12 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  gamma2p <- c(gamma2[1], gamma2[2] + epsilon)
  gamma2m <- c(gamma2[1], gamma2[2] - epsilon)
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2p, beta1, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2m, beta1, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2p
                          , beta2, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2m
                          , beta2, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialB01 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  beta1p <- c(beta1[1] + epsilon, beta1[2])
  beta1m <- c(beta1[1] - epsilon, beta1[2])
  delta1p <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1p, prec1, tau, p, x)))
  delta1m <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1m, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1p)), 1, function(a) SolveDelta2(gamma2, beta1p, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1m)), 1, function(a) SolveDelta2(gamma2, beta1m, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1p, beta1p, prec1, delta2p, beta2, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1m, beta1m, prec1, delta2m, beta2, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialB11 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  beta1p <- c(beta1[1], beta1[2] + epsilon)
  beta1m <- c(beta1[1], beta1[2] - epsilon)
  delta1p <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1p, prec1, tau, p, x)))
  delta1m <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1m, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1p)), 1, function(a) SolveDelta2(gamma2, beta1p, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1m)), 1, function(a) SolveDelta2(gamma2, beta1m, beta2, prec1, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1p, beta1p, prec1, delta2p, beta2, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1m, beta1m, prec1, delta2m, beta2, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialB02 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  beta2p <- c(beta2[1] + epsilon, beta2[2], beta2[3])
  beta2m <- c(beta2[1] - epsilon, beta2[2], beta2[3])
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2p, prec1, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2m, prec1, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2p, beta2p, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2m, beta2m, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialB12 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  beta2p <- c(beta2[1], beta2[2] + epsilon, beta2[3])
  beta2m <- c(beta2[1], beta2[2] - epsilon, beta2[3])
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2p, prec1, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2m, prec1, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2p, beta2p, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2m, beta2m, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialB22 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  beta2p <- c(beta2[1], beta2[2], beta2[3] + epsilon)
  beta2m <- c(beta2[1], beta2[2], beta2[3] - epsilon)
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2p, prec1, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2m, prec1, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2p, beta2p, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2m, beta2m, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialP11 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  prec1p <- c(prec1[1] + epsilon, prec1[2])
  prec1m <- c(prec1[1] - epsilon, prec1[2])
  delta1p <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1p, tau, p, x)))
  delta1m <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1m, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1p)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1p, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1m)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1m, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1p, beta1, prec1p, delta2p, beta2, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1m, beta1, prec1m, delta2m, beta2, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialP01 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  prec1p <- c(prec1[1], prec1[2] + epsilon)
  prec1m <- c(prec1[1], prec1[2] - epsilon)
  delta1p <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1p, tau, p, x)))
  delta1m <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1m, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1p)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1p, prec2, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1m)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1m, prec2, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1p, beta1, prec1p, delta2p, beta2, prec2, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1m, beta1, prec1m, delta2m, beta2, prec2, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialP12 <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  prec2p <- prec2 + epsilon
  prec2m <- prec2 - epsilon
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2p, lambda, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2m, lambda, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2p, beta2, prec2p, lambda)
  llm <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2m, beta2, prec2m, lambda)
  return((llp - llm) / 2/ epsilon)
}

PartialLambda <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  lambdap <- lambda + epsilon
  lambdam <- lambda - epsilon
  delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, p, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2, lambdap, tau, p, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2, lambdam, tau, p, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2p, beta2, prec2, lambdap)
  llm <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2m, beta2, prec2, lambdam)
  return((llp - llm) / 2/ epsilon)
}

PartialPhi <- function(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, p, tau, y, R, x){
  epsilon <- 0.003
  n <- length(x)
  phi1 <- p + epsilon
  phi2 <- p - epsilon
  delta1p <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, phi1, x)))
  delta1m <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, phi2, x)))
  delta2p <- apply(as.matrix(cbind(x, delta1p)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2, lambda, tau, phi1, a[1], a[2]))
  delta2m <- apply(as.matrix(cbind(x, delta1m)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2, lambda, tau, phi2, a[1], a[2]))
  llp <- NegLogLikelihood(y, R, x, delta1p, beta1, prec1, delta2p, beta2, prec2, lambda) - sum(R) * log(phi1) - (n - sum(R)) * log(1 - phi1)
  llm <- NegLogLikelihood(y, R, x, delta1m, beta1, prec1, delta2m, beta2, prec2, lambda) - sum(R) * log(phi2) - (n - sum(R)) * log(1 - phi2)
  return((llp - llm) / 2/ epsilon)
}


###############
## Main Function 
###############

BiQRGradient <- function(y, R, x, tau, niter = 1000){
  phi <- 0.5
  n <- length(R)
  gamma1 <- c(0, 0)
  gamma2 <- c(0, 0)
  beta1 <- c(0, 0) ## beta_1^(1) 
  beta2 <- c(0, 0, 0) ## beta_2^(0) == 0 for MAR
  prec1 <- c(1, 1)
  prec2 <- 1
  lambda <- 1 ## lambda == 1 for MAR
  delta1 <- rep(0, n)
  delta2 <- rep(0, n)
  ll0 <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2, beta2, prec2, lambda)

  paramsave <- c(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, ll0)
  
  dif <- 1
  alpha <- 0.0003
  iter <- 0

  ## save

  ## progress bar
  pb <- txtProgressBar(min = 0, max = niter, style = 3)

  ## BEGIN GRADIENT DESCENT
  while (dif > 10^-5 & iter < niter) {
    ## GET PARTIAL
    pg01 <- PartialG01(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    pg11 <- PartialG11(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    pg02 <- PartialG02(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    pg12 <- PartialG12(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    pb01 <- PartialB01(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    pb11 <- PartialB11(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    ## pb02 <- PartialB02(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    ## pb12 <- PartialB12(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    ## pb22 <- PartialB22(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    pp11 <- PartialP11(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    pp01 <- PartialP01(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    pp12 <- PartialP12(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
##    plambda <- PartialLambda(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)
    pphi <- PartialPhi(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, phi, tau, y, R, x)

    ## MOVE
    gamma1[1] <- gamma1[1] - alpha * pg01
    gamma1[2] <- gamma1[2] - alpha * pg11
    gamma2[1] <- gamma2[1] - alpha * pg02
    gamma2[2] <- gamma2[2] - alpha * pg12
    beta1[1] <- beta1[1] - alpha * pb01
    beta1[2] <- beta1[2] - alpha * pb11
    ## beta2[1] <- beta2[1] - alpha * pb02
    ## beta2[2] <- beta2[2] - alpha * pb12
    ## beta2[3] <- beta2[3] - alpha * pb22
##    prec1[1] <- max(prec1[1] - alpha * pp01, 0.01)
    prec1[1] <- max(prec1[1] - 0.00003 * pp01, 0.01)
##    prec1[2] <- max(prec1[2] - alpha * pp11, 0.01)
    prec1[2] <- max(prec1[2] - 0.00003 * pp11, 0.01)
##    prec2 <- max(prec2 - alpha * pp12, 0.01)
    prec2 <- max(prec2 - 0.00003 * pp12, 0.01)
    ## lambda <- lambda - alpha * plambda
    phi <- max(min(phi - alpha * pphi, 0.99), 0.01)

    ## CALCULATE DIF
    delta1 <- unlist(sapply(as.list(x), function(x) SolveDelta1(gamma1, beta1, prec1, tau, phi, x)))
    delta2 <- apply(as.matrix(cbind(x, delta1)), 1, function(a) SolveDelta2(gamma2, beta1, beta2, prec1, prec2, lambda, tau, phi, a[1], a[2]))
    ll <- NegLogLikelihood(y, R, x, delta1, beta1, prec1, delta2, beta2, prec2, lambda) - sum(R) * log(phi) - (n - sum(R)) * log(1 - phi)
    dif <- abs(ll - ll0)
    ll0 <- ll

    ## SAVE
    paramsave <- rbind(paramsave, c(gamma1, gamma2, beta1, beta2, prec1, prec2, lambda, ll))
    
    iter <- iter + 1

    ## PROGRESSBAR
    setTxtProgressBar(pb, iter)
  }

  param <- list(gamma1 = gamma1, gamma2 = gamma2, beta1 = beta1, beta2 = beta2, prec1 = prec1, prec2 = prec2, lambda = lambda, phi = phi, iter = iter, tau = tau, paramsave = paramsave)

  return(param)
  
}

myplot <- function(mod){
  plot(y[,1] ~x)
  abline(mod$gamma1)
  plot(y[,2] ~x)
  abline(mod$gamma2)
  apply(mod$paramsave, 2, function(x) plot(ts(x)))
}

