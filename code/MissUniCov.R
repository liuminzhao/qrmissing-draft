# univariate
# missing
# general
# with covariates
# Time-stamp: <liuminzhao 03/16/2013 11:04:13>

# simulate data
# y | S = 1 ~ N(\delta + b0 + x*b1, sigma1)
# y | S = 0 ~ N(\delta - b0 - x*b1, sigma2)
# P(S = 1) = pi

# y | S = 1 ~ N(b01 + x*b11, sigma1)
# y | S = 0 ~ N(b00 + x*b10, sigma2)

library(rootSolve)

n <- 200
p <- 0.5
S <- rbinom(n, 1, p)
b01 <- 1
b00 <- -1
b11 <- 2
b10 <- -2
sigma1 <- 1
sigma0 <- 1
x <- rnorm(n)
y <- rep(0, n)
for (i in 1:n){
  if (S[i] == 1) {
    y[i] <- rnorm(1, b01 + x[i] * b11, sigma1)
  } else {
    y[i] <- rnorm(1, b00 + x[i] * b10, sigma0)
  }
}
y1 <- y[S == 1]
y0 <- y[S == 0]

# mean regression of y on x
summary(lm(y ~ x))
plot(y ~ x)

tau <- 0.9
phi <- 0.5

## initial
beta <- c(1, 2) # beta^(1)
gamma <- c(0, 0)
sigma <- c(sigma1, sigma0)
delta <- rep(0, n)
ll <- LogLikelihood(y, S, x, delta, beta, sigma)

## MCMC paramters
nsave <- 10000
nburn <- nsave/2
niter <- nburn + nsave

## save parameter
betasave <- matrix(0, 2, nsave)
gammasave <- matrix(0, 2, nsave)
sigmasave <- matrix(0, 2, nsave)

## TUNING
att <- 0
acc <- 0
tune <- 1
thres <- 0.245

## begin MCMC
current <- proc.time()[3]
pb <- txtProgressBar(min = 0, max = niter, style = 3)
for (iter in 1:niter) {
  ## update gamma
  att <- att + 1
  gammac <- gamma + rnorm(2, sd = tune)
  deltac <- unlist(sapply(as.list(x), function(x) SolveDelta(gammac, beta, sigma, tau, phi, x)))
  llc <- LogLikelihood(y, S, x, deltac, beta, sigma)
  ratio <- llc - ll
  if (runif(1) < exp(ratio)) {
    gamma <- gammac
    ll <- llc
    acc <- acc + 1
  }
  ## TUNING
  if (iter <= nburn){
    if (att >= 100){
      rate <- acc / att
      if (rate > thres){
        tune <- tune * 1.1
      } else {
        tune <- tune / 1.1
      }
      tune <- min(tune, 100)
      tune <- max(tune, 0.01)
      att <- acc <- 0
    }
  }
  ## save
  if (iter > nburn) {
    gammasave[, iter - nburn] <- gamma
  }
  ## Update Pb
  setTxtProgressBar(pb, iter)
}
close(pb)
cat('Time cost ' , proc.time()[3] - current)


###############
## diagnostic 
###############

plot(ts(gammasave[1, ]))
plot(ts(gammasave[2, ]))

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param gamma 
##' @param beta 
##' @param sigma 
##' @param tau 
##' @param p 
##' @param y 
##' @param x 
##' @return target function
##' @author Minzhao Liu
TargetEqn <- function(delta, gamma, beta, sigma, tau, p, x){
  quan <- gamma[1] + gamma[2] * x
  lp <- beta[1] + beta[2] * x
  return(tau - p * pnorm((quan - delta - lp)/sigma[1]) - (1 - p) * pnorm((quan - delta + lp)/sigma[2]))
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param gamma 
##' @param beta 
##' @param sigma 
##' @param tau 
##' @param p 
##' @param y 
##' @param x 
##' @return delta
##' @author Minzhao Liu
SolveDelta <- function(gamma, beta, sigma, tau, p, x){
  return(uniroot.all(TargetEqn, c(-30, 30), tol = 0.0001, gamma = gamma, beta = beta, sigma = sigma, tau = tau, p = p, x = x))
}

LogLikelihood <- function(y, S, x, delta, beta, sigma){
  ll <- 0
  mu1 <- delta + beta[1] + beta[2] * x
  mu0 <- delta - beta[1] - beta[2] * x
  ll1 <- dnorm(y, mean = mu1, sd = sigma1, log = T)
  ll0 <- dnorm(y, mean = mu0, sd = sigma0, log = T)
  ll <- c(ll1[S == 1], ll0[S == 0])
  return(sum(ll))
}



Test <- function(){
  ## test TargetEqn and SolveDelta
  
  ## beta0 == beta1 = 0 and sigma1 = sigma0
  ## delta = x \gamma - sigma*qnorm(tau)
  gamma <- c(1, 2)
  s <- 3
  sigma <- c(s, s)
  tau <- 0.5
  x <- 5
  if (isTRUE(all.equal(SolveDelta(gamma, c(0,0), sigma, tau, 0.9, x), gamma[1] + gamma[2] * x - s * qnorm(tau), tol = 0.00001))) print(TRUE)

  gamma <- c(1, -2)
  s <- 2
  sigma <- c(s, s)
  tau <- 0.9
  x <- 3
  if (isTRUE(all.equal(SolveDelta(gamma, c(0,0), sigma, tau, 0.5, x), gamma[1] + gamma[2] * x - s * qnorm(tau), tol = 0.00001))) print(TRUE)

}
