#!/bin/Rscript
##' 2013/06/27 simulation on Laplace distribution
##' 2013/07/04 new simulation
##' 2013/07/10 try small value and add Bottai's algorithm

sink('sim-lp-0710-ex.txt')
rm(list = ls())
source('BiMLESigma.R')
source('sendEmail.R')
source('Bottai.R')
library(quantreg)
library(xtable)
library(doMC)
registerDoMC()
options(cores = 8)
set.seed(1)

###############
## PARAMETER
###############
n <- 500
p <- 0.5
alpha <- 0.5
tau <- 1

###############
## SIMULATION
###############
boot <- 1000

start <- proc.time()[3]

result <- foreach(icount(boot), .combine = rbind) %dopar% {
  R <- rbinom(n, 1, p)
  x <- runif(n, 0, 2)
  y <- matrix(0, n, 2)
  w <- rgamma(n, 1, rate = tau)
  for (i in 1:n){
    if (R[i] == 1){
      y[i, 1] <- 2 + x[i] + (1 + alpha*x[i])*rnorm(1, 0, sd = sqrt(8*w[i]/tau))
      y[i, 2] <- 1 - x[i] - y[i, 1]*1/2 + rnorm(1,sd = sqrt(6*w[i]/tau))*(1 + alpha*x[i])
    } else {
      y[i, 1] <- -2 - x[i] + (1 + alpha*x[i])* rnorm(1, sd = sqrt(8*w[i]/tau))
      y[i, 2] <- 1 - x[i] - y[i, 1]/2 + rnorm(1, sd = sqrt(6*w[i]/tau))*(1 + alpha*x[i])
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

  mod1mm <- coef(mod1)
  mod3mm <- coef(mod3)
  mod5mm <- coef(mod5)
  mod7mm <- coef(mod7)
  mod9mm <- coef(mod9)

  mod1rq <- as.vector(rq(y[,1]~x, tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)
  mod2rq <- as.vector(rq(y[,2][R==1]~x[R==1], tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)

  mod1b <- Bottai(y, R, X, tau = 0.1)
  mod3b <- Bottai(y, R, X, tau = 0.3)
  mod5b <- Bottai(y, R, X, tau = 0.5)
  mod7b <- Bottai(y, R, X, tau = 0.7)
  mod9b <- Bottai(y, R, X, tau = 0.9)

  ans <- c(mod1mm[1,], mod3mm[1,], mod5mm[1,], mod7mm[1,], mod9mm[1,],
           mod1mm[2,], mod3mm[2,], mod5mm[2,], mod7mm[2,], mod9mm[2,],
           mod1rq, mod2rq,
           mod1b[,1], mod3b[,1], mod5b[,1], mod7b[,1], mod9b[,1],
           mod1b[,2], mod3b[,2], mod5b[,2], mod7b[,2], mod9b[,2])
}

write.table(result, file = "sim-lp-error-0710.txt", row.names = F, col.names = F)
sendEmail(subject="simulation-lp-MAR", text="done", address="liuminzhao@gmail.com")

###############
## AFTER SIM; GET SUMMARY
###############

###############
## TRUE VALUE
###############

pLD <- function(x, tau){
  return(ifelse(x < 0, exp(tau * x)/2, 1 - exp(-tau * x)/2))
}

quan1 <- function(y, x, tau, quan){
  return(quan - .5*pLD(y - 2- x, tau = tau/(1 + alpha*x) ) - .5*pLD(y + 2 +x, tau = tau/(1 + alpha*x)))
}

SolveQuan1 <- function(x, tau, quan){
  uniroot(quan1, c(-30, 30), x = x, tau = tau, quan = quan)$root
}

quan2 <- function(y, x, tau, quan){
  return(quan - .5*pLD(y + (1.5*x), tau = tau/(1 + alpha*x)) - .5*pLD(y - (2 - 0.5*x),tau = tau/(1 + alpha*x)))
}

SolveQuan2 <- function(x, tau, quan){
  uniroot(quan2, c(-30, 30), x = x, tau = tau, quan = quan)$root
}

xsim <- seq(0, 2, len = 100)
y15 <- sapply(xsim, function(x) SolveQuan1(x, tau, 0.5))
y19 <- sapply(xsim, function(x) SolveQuan1(x, tau, 0.9))
y17 <- sapply(xsim, function(x) SolveQuan1(x, tau, 0.7))
y13 <- sapply(xsim, function(x) SolveQuan1(x, tau, 0.3))
y11 <- sapply(xsim, function(x) SolveQuan1(x, tau, 0.1))

q11 <- lm(y11~xsim)$coef
q13 <- lm(y13~xsim)$coef
q15 <- lm(y15~xsim)$coef
q17 <- lm(y17~xsim)$coef
q19 <- lm(y19~xsim)$coef

y25 <- sapply(xsim, function(x) SolveQuan2(x, tau, 0.5))
y29 <- sapply(xsim, function(x) SolveQuan2(x, tau, 0.9))
y27 <- sapply(xsim, function(x) SolveQuan2(x, tau, 0.7))
y23 <- sapply(xsim, function(x) SolveQuan2(x, tau, 0.3))
y21 <- sapply(xsim, function(x) SolveQuan2(x, tau, 0.1))

q21 <- lm(y21~xsim)$coef
q23 <- lm(y23~xsim)$coef
q25 <- lm(y25~xsim)$coef
q27 <- lm(y27~xsim)$coef
q29 <- lm(y29~xsim)$coef

result <- read.table('sim-lp-error-0710.txt')
trueq <- c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
trueq <- rep(trueq, 3)
mse <- rep(0, 60)
for (i in 1:60){
  mse[i] <- mean((result[,i] - trueq[i])^2, trim = 0.05)
}

mseh2 <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))
mserq <- rbind(matrix(mse[21:30], 2, 5), matrix(mse[31:40], 2, 5))
msebb <- rbind(matrix(mse[41:50], 2, 5), matrix(mse[51:60], 2, 5))
mytbl <- cbind(mseh2[,1], mserq[,1], msebb[,1], mseh2[,2], mserq[,2], msebb[,2], mseh2[,3], mserq[,3], msebb[,3], mseh2[,4], mserq[,4], msebb[,4], mseh2[,5], mserq[,5], msebb[,5])

colnames(mytbl) <- rep(c('MM', 'RQ', 'Bottai'), 5)

print(xtable(mytbl))
print(mytbl)
cat("Time: ", proc.time()[3] - start, '\n')
sink()
