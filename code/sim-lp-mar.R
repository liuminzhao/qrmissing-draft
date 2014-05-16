#!/bin/Rscript
##' 2013/06/27 simulation on Laplace distribution
##' 2013/07/04 new simulation
##' 2013/07/10 try small value and add Bottai's algorithm
##' 2013/07/13 try 4 +- x
##' 2013/07/16 try +- 4 +- 2x
##' 2013/08/01 on QRMissingBi
##' 2013/08/01 test on QRMissingBi.R

sink('sim-lp-mar-0516.txt')
rm(list = ls())
library(xtable)
library(qrmissing)
library(doMC)
source('sendEmail.R')
source('Bottai.R')
registerDoMC()
options(cores = 10)
set.seed(1)

###############
## PARAMETER
###############
n <- 200
p <- 0.5
alpha <- 0
tau <- 1

###############
## SIMULATION
###############
boot <- 100

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

  mod1 <- QRMissingBiMixMLE(y, R, X, tau = 0.1, K = 2, model = 'slope')
  mod3 <- QRMissingBiMixMLE(y, R, X, tau = 0.3, K = 2, model = 'slope')
  mod5 <- QRMissingBiMixMLE(y, R, X, tau = 0.5, K = 2, model = 'slope')
  mod7 <- QRMissingBiMixMLE(y, R, X, tau = 0.7, K = 2, model = 'slope')
  mod9 <- QRMissingBiMixMLE(y, R, X, tau = 0.9, K = 2, model = 'slope')

  mod1mm <- rbind(coef(mod1)$gamma1, coef(mod1)$gamma2)
  mod3mm <- rbind(coef(mod3)$gamma1, coef(mod3)$gamma2)
  mod5mm <- rbind(coef(mod5)$gamma1, coef(mod5)$gamma2)
  mod7mm <- rbind(coef(mod7)$gamma1, coef(mod7)$gamma2)
  mod9mm <- rbind(coef(mod9)$gamma1, coef(mod9)$gamma2)

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

write.table(result, file = "sim-lp-mar-0516-result.txt", row.names = F, col.names = F)
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
  return(quan - .5*pLD(y - 2 - x, tau = tau/(1 + alpha*x) ) - .5*pLD(y + 2 + x, tau = tau/(1 + alpha*x)))
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


result <- read.table('sim-lp-mar-0516-result.txt')
trueq <- c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
trueq <- rep(trueq, 3)

## MSE
mse <- rep(0, 60)
for (i in 1:60){
  mse[i] <- mean((result[,i] - trueq[i])^2, trim = 0.05)
}
mseh2 <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))
mserq <- rbind(matrix(mse[21:30], 2, 5), matrix(mse[31:40], 2, 5))
msebb <- rbind(matrix(mse[41:50], 2, 5), matrix(mse[51:60], 2, 5))
mytbl <- cbind(mseh2[,1], mserq[,1], msebb[,1], mseh2[,2], mserq[,2], msebb[,2], mseh2[,3], mserq[,3], msebb[,3], mseh2[,4], mserq[,4], msebb[,4], mseh2[,5], mserq[,5], msebb[,5])
colnames(mytbl) <- rep(c('MM', 'RQ', 'BB'), 5)
print(xtable(mytbl))
print(mytbl)

## BIAS
bias <- rep(0, 60)
for (i in 1:60){
  bias[i] <- mean((result[,i] - trueq[i]), trim = 0.05)
}
biash2 <- rbind(matrix(bias[1:10], 2, 5), matrix(bias[11:20], 2, 5))
biasrq <- rbind(matrix(bias[21:30], 2, 5), matrix(bias[31:40], 2, 5))
biasbb <- rbind(matrix(bias[41:50], 2, 5), matrix(bias[51:60], 2, 5))
mytbl2 <- cbind(biash2[,1], biasrq[,1], biasbb[,1], biash2[,2], biasrq[,2], biasbb[,2], biash2[,3], biasrq[,3], biasbb[,3], biash2[,4], biasrq[,4], biasbb[,4], biash2[,5], biasrq[,5], biasbb[,5])
colnames(mytbl2) <- rep(c('MM', 'RQ', 'BB'), 5)
print(xtable(mytbl2))
print(mytbl2)

## Efficiency
efficiency <- apply(result[,1:60], 2, var)
efficiencyh2 <- rbind(matrix(efficiency[1:10], 2, 5), matrix(efficiency[11:20], 2, 5))
efficiencyrq <- rbind(matrix(efficiency[21:30], 2, 5), matrix(efficiency[31:40], 2, 5))
efficiencybb <- rbind(matrix(efficiency[41:50], 2, 5), matrix(efficiency[51:60], 2, 5))
mytbl3 <- cbind(efficiencyrq/efficiencyh2, efficiencybb/efficiencyh2)
colnames(mytbl3) <- rep(c('RQ', 'BB'), 5)
print(xtable(mytbl3))

## MCSE
mcse <- rep(0, 60)
for (i in 1:60){
  mcse[i] <- sd((result[,i] - trueq[i])^2)/sqrt(boot)
}
mcseh2 <- rbind(matrix(mcse[1:10], 2, 5), matrix(mcse[11:20], 2, 5))
mcserq <- rbind(matrix(mcse[21:30], 2, 5), matrix(mcse[31:40], 2, 5))
mcsebb <- rbind(matrix(mcse[41:50], 2, 5), matrix(mcse[51:60], 2, 5))
mytbl4 <- cbind(mcseh2[,1], mcserq[,1], mcsebb[,1], mcseh2[,2], mcserq[,2], mcsebb[,2], mcseh2[,3], mcserq[,3], mcsebb[,3], mcseh2[,4], mcserq[,4], mcsebb[,4], mcseh2[,5], mcserq[,5], mcsebb[,5])
colnames(mytbl4) <- rep(c('MM', 'RQ', 'BB'), 5)
print(xtable(mytbl4))
print(mytbl4)

## combine mse and MCSE
msemcse <- matrix(0, 4, 30)
for (i in 1:15){
  msemcse[, i*2 - 1] <- mytbl[, i]
  msemcse[, i*2] <- mytbl4[, i]
}
print(msemcse)
print(xtable(msemcse))

cat("Time: ", proc.time()[3] - start, '\n')
sink()
