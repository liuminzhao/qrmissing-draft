##' Time-stamp: <liuminzhao 06/24/2013 16:08:28>
##' Simulation for paper,
##' T error
##' 2013/06/24
##' MAR

rm(list = ls())
source('BiMLESigma.R')
source('sendEmail.R')
library(quantreg)
library(xtable)
library(doMC)
registerDoMC()
options(cores = 8)
set.seed(1)


###############
## TRUE VALUE
###############

quan1 <- function(y, x, tau){
  return(tau - .5*pt(y - 1- x, df = 3) - .5*pt(y+1+x,df = 3))
}

SolveQuan1 <- function(x, tau){
  uniroot(quan1, c(-30, 30), x = x, tau = tau)$root
}

x <- seq(0, 2, len = 100)
y15 <- sapply(x, function(x) SolveQuan1(x, 0.5))
y19 <- sapply(x, function(x) SolveQuan1(x, 0.9))
y17 <- sapply(x, function(x) SolveQuan1(x, 0.7))
y13 <- sapply(x, function(x) SolveQuan1(x, 0.3))
y11 <- sapply(x, function(x) SolveQuan1(x, 0.1))

q11 <- lm(y11~x)$coef
q13 <- lm(y13~x)$coef
q15 <- lm(y15~x)$coef
q17 <- lm(y17~x)$coef
q19 <- lm(y19~x)$coef

quan2 <- function(y, x, tau){
  return(tau - .5*pt(y - 1.5- 1.5*x, df = 3) - .5*pt(y - 0.5 - 0.5*x,df = 3))
}

SolveQuan2 <- function(x, tau){
  uniroot(quan2, c(-30, 30), x = x, tau = tau)$root
}

x <- seq(0, 2, len = 100)
y25 <- sapply(x, function(x) SolveQuan2(x, 0.5))
y29 <- sapply(x, function(x) SolveQuan2(x, 0.9))
y27 <- sapply(x, function(x) SolveQuan2(x, 0.7))
y23 <- sapply(x, function(x) SolveQuan2(x, 0.3))
y21 <- sapply(x, function(x) SolveQuan2(x, 0.1))

q21 <- lm(y21~x)$coef
q23 <- lm(y23~x)$coef
q25 <- lm(y25~x)$coef
q27 <- lm(y27~x)$coef
q29 <- lm(y29~x)$coef

###############
## PARAMETER
###############
n <- 500
p <- 0.5

###############
## SIMULATION
###############
boot <- 8

result <- foreach(icount(boot), .combine = rbind) %dopar% {
  R <- rbinom(n, 1, p)
  x <- runif(n, 0, 2)
  y <- matrix(0, n, 2)

  v <- 3
  w <- rgamma(n, v/2, rate = v/2)
  winv <- 1/w


  for (i in 1:n){
    if (R[i] == 1){
      y[i, 1] <- rnorm(1, 1 + x[i], sd = sqrt(winv[i]))
      y[i, 2] <- rnorm(1, 1 + x[i] + y[i, 1]*1/2, sd = sqrt(3*winv[i]/4))
    } else {
      y[i, 1] <- rnorm(1, -1 - x[i], sd = sqrt(winv[i]))
      y[i, 2] <- rnorm(1, 1 + x[i] + y[i, 1]/2, sd = sqrt(3*winv[i]/4))
    }
  }

  X <- matrix(0, n, 2)
  X[,1] <- 1
  X[,2] <- x

  mod1 <- BiQRGradient(y, R, X, tau = 0.1, method = 'heter2')
  mod3 <- BiQRGradient(y, R, X, tau = 0.3, method = 'heter2')
  mod5 <- BiQRGradient(y, R, X, tau = 0.5, method = 'heter2')
  mod7 <- BiQRGradient(y, R, X, tau = 0.7, method = 'heter2')
  mod9 <- BiQRGradient(y, R, X, tau = 0.9, method = 'heter2')

  mod1mm <- mysummary(mod1)
  mod3mm <- mysummary(mod3)
  mod5mm <- mysummary(mod5)
  mod7mm <- mysummary(mod7)
  mod9mm <- mysummary(mod9)

  mod1rq <- as.vector(rq(y[,1]~x, tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)
  mod2rq <- as.vector(rq(y[,2][R==1]~x[R==1], tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)

  ans <- c(mod1mm[1,], mod3mm[1,], mod5mm[1,], mod7mm[1,], mod9mm[1,],
           mod1mm[2,], mod3mm[2,], mod5mm[2,], mod7mm[2,], mod9mm[2,], mod1rq, mod2rq)
}

write.table(result, file = "sim-t-error.txt", row.names = F, col.names = F)
sendEmail(subject="simulation-t-MAR", text="done", address="liuminzhao@gmail.com")
