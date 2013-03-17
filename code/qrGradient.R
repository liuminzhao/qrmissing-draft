## univariate
## missing
## general
## with covariates
## MLE with gradient descent method
## Time-stamp: <liuminzhao 03/16/2013 23:45:56>

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
  ll1 <- dnorm(y, mean = mu1, sd = sigma1, log = T)
  ll0 <- dnorm(y, mean = mu0, sd = sigma0, log = T)
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


###############
## Test 
###############

n <- 200
p <- 0.5
S <- rbinom(n, 1, p)
b01 <- 1
b00 <- -1
b11 <- 2
b10 <- -2
sigma1 <- 1
sigma0 <- 1
x <- runif(n, 0, 2)
## x <- rnorm(n)
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

## mean regression of y on x
summary(lm(y ~ x))
plot(y ~ x)

tau <- 0.1
phi <- 0.5
beta <- c(1, 2) # beta^(1)
gamma <- c(0, 0)
sigma <- c(sigma1, sigma0)
delta <- rep(0, n)
ll0 <- LogLikelihood(y, S, x, delta, beta, sigma)
dif <- 1
alpha <- 0.001
iter <- 0
llsave <- ll0


###############
## BEGIN GRADIENT DESCENT 
###############
while (dif > 10^-3 & iter < 1000) {
  pg0 <- PartialG0(gamma, beta)
  pg1 <- PartialG1(gamma, beta)
  pb0 <- PartialB0(gamma, beta)
  pb1 <- PartialB1(gamma, beta)
  gamma[1] <- gamma[1] - alpha * pg0
  gamma[2] <- gamma[2] - alpha * pg1
  beta[1] <- beta[1] - alpha * pb0
  beta[2] <- beta[2] - alpha * pb1
  delta <- unlist(sapply(as.list(x), function(x) SolveDelta(gamma, beta, sigma, tau, phi, x)))
  ll <- LogLikelihood(y, S, x, delta, beta, sigma)
  dif <- abs(ll - ll0)
  ll0 <- ll
  llsave <- append(llsave, ll)
  iter <- iter + 1
}
cat(gamma, beta, sigma, '\n')
plot(ts(llsave))


###############
## save png 
###############

png('../image/mle1.png')
plot(y ~x)
abline(c(1.93, 2.038), col = 2)
abline(c(0.91, 1.977), col = 3)
abline(c(0.06, -0.012), col = 4)
abline(c(-0.733, -2.06), col = 5)
abline(c(-1.79, -2.081), col = 6)
legend('topleft', c('0.1', '0.3', '0.5', '0.7', '0.9'), col = 6:2, lty = rep(1, 5))
dev.off()
