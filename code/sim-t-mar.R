#!/bin/Rscript
##' Time-stamp: <liuminzhao 01/10/2015 22:46:46>
##' Simulation for paper,
##' T error
##' 2013/06/24
##' MAR
##' 2013/07/02 new t3 error
##' 2013/07/07 smaller version of t error like normal
##' 2013/07/11 keep alpha = 0.5 instead of 0.3 and do MM, RQ, BB together
##' 2013/08/01 test on QRMissingBi.R

sink('sim-t-mar-0101.txt')
rm(list = ls())
library(xtable)
library(qrmissing)
## library(doMC)
source('sendEmail.R')
source('Bottai.R')
source('QR_Panel.R') # QR with longitudinal
source('Lipsitz.R') # Lipsitz
## registerDoMC()
options(cores = 1)
set.seed(1)

###############
## PARAMETER
###############
n <- 200
p <- 0.5
alpha <- 0

###############
## SIMULATION
###############
boot <- 20

start <- proc.time()[3]

## result <- foreach(icount(boot), .combine = rbind) %dopar% {
for (icount in 1:boot){
  R <- rbinom(n, 1, p)
  x <- runif(n, 0, 2)
  y <- matrix(0, n, 2)

  v <- 3
  w <- rgamma(n, v/2, rate = v/2)
  winv <- 1/w

  for (i in 1:n){
    if (R[i] == 1){
    y[i, 1] <-  2 + x[i] + (1 + alpha*x[i])*rnorm(1, 0, sd = sqrt(winv[i]))
    y[i, 2] <-  1 - x[i] - y[i, 1]*1/2  +  (1 + alpha*x[i])*rnorm(1, 0, sd = sqrt(3*winv[i]/4))
  } else {
    y[i, 1] <-  -2 - x[i] +(1 + alpha*x[i])*rnorm(1, 0, sd = sqrt(winv[i]))
    y[i, 2] <-  1 - x[i] - y[i, 1]*1/2  +  (1 + alpha*x[i])*rnorm(1, 0, sd = sqrt(3*winv[i]/4))
    }
  }

  X <- matrix(0, n, 2)
  X[,1] <- 1
  X[,2] <- x

  ## M1

  mod1 <- QRMissingBi(y ~ x, R, tau = 0.1, model = 'slope')
  mod3 <- QRMissingBi(y ~ x, R, tau = 0.3, model = 'slope')
  mod5 <- QRMissingBi(y ~ x, R, tau = 0.5, model = 'slope')
  mod7 <- QRMissingBi(y ~ x, R, tau = 0.7, model = 'slope')
  mod9 <- QRMissingBi(y ~ x, R, tau = 0.9, model = 'slope')

  mod1mm <- rbind(coef(mod1)$gamma1, coef(mod1)$gamma2)
  mod3mm <- rbind(coef(mod3)$gamma1, coef(mod3)$gamma2)
  mod5mm <- rbind(coef(mod5)$gamma1, coef(mod5)$gamma2)
  mod7mm <- rbind(coef(mod7)$gamma1, coef(mod7)$gamma2)
  mod9mm <- rbind(coef(mod9)$gamma1, coef(mod9)$gamma2)

  BIC1 <- mod1$BIC
  BIC3 <- mod3$BIC
  BIC5 <- mod5$BIC
  BIC7 <- mod7$BIC
  BIC9 <- mod9$BIC

  ## M2 BIC choose

  for (k in 2:2) {
      mod1tmp <- QRMissingBiMixMLE(y~ x,R,tau = 0.1, K = k, model = 'slope')
      mod3tmp <- QRMissingBiMixMLE(y~ x,R,tau = 0.3, K = k, model = 'slope')
      mod5tmp <- QRMissingBiMixMLE(y~ x,R,tau = 0.5, K = k, model = 'slope')
      mod7tmp <- QRMissingBiMixMLE(y~ x,R,tau = 0.7, K = k, model = 'slope')
      mod9tmp <- QRMissingBiMixMLE(y~ x,R,tau = 0.9, K = k, model = 'slope')
      BIC1tmp <- mod1tmp$BIC
      BIC3tmp <- mod3tmp$BIC
      BIC5tmp <- mod5tmp$BIC
      BIC7tmp <- mod7tmp$BIC
      BIC9tmp <- mod9tmp$BIC
      if (BIC1tmp < BIC1) {
          BIC1 <- BIC1tmp
          mod1 <- mod1tmp
      }
      if (BIC3tmp < BIC3) {
          BIC3 <- BIC3tmp
          mod3 <- mod3tmp
      }
      if (BIC5tmp < BIC5) {
          BIC5 <- BIC5tmp
          mod5 <- mod5tmp
      }
      if (BIC7tmp < BIC7) {
          BIC7 <- BIC7tmp
          mod7 <- mod7tmp
      }
      if (BIC9tmp < BIC9) {
          BIC9 <- BIC9tmp
          mod9 <- mod9tmp
      }
  }

  mod1m2coef <- rbind(coef(mod1)$gamma1, coef(mod1)$gamma2)
  mod3m2coef <- rbind(coef(mod3)$gamma1, coef(mod3)$gamma2)
  mod5m2coef <- rbind(coef(mod5)$gamma1, coef(mod5)$gamma2)
  mod7m2coef <- rbind(coef(mod7)$gamma1, coef(mod7)$gamma2)
  mod9m2coef <- rbind(coef(mod9)$gamma1, coef(mod9)$gamma2)

  ## RQ and BZ
  modQR <- QR.Panel(X, y, R, w = rep(1/5, 5), taus = c(0.1, 0.3, 0.5, 0.7, 0.9))

  ## Bottai
  mod1b <- Bottai(y, R, X, tau = 0.1)
  mod3b <- Bottai(y, R, X, tau = 0.3)
  mod5b <- Bottai(y, R, X, tau = 0.5)
  mod7b <- Bottai(y, R, X, tau = 0.7)
  mod9b <- Bottai(y, R, X, tau = 0.9)
  ## Lipsitz
  modLip <- Lipsitz(X, y, R, tau = c(0.1, 0.3, 0.5, 0.7, 0.9))


  ## All results together

  ans <- c(mod1mm[1,], mod3mm[1,], mod5mm[1,], mod7mm[1,], mod9mm[1,],
           mod1mm[2,], mod3mm[2,], mod5mm[2,], mod7mm[2,], mod9mm[2,],
           mod1m2coef[1,], mod3m2coef[1,], mod5m2coef[1,], mod7m2coef[1,], mod9m2coef[1,],
           mod1m2coef[2,], mod3m2coef[2,], mod5m2coef[2,], mod7m2coef[2,], mod9m2coef[2,],
           c(modQR$coef[1:2, ]), c(modQR$coef[3:4, ]),
           c(modLip$coef[1:2, ]), c(modLip$coef[3:4, ]),
           mod1b[,1], mod3b[,1], mod5b[,1], mod7b[,1], mod9b[,1],
           mod1b[,2], mod3b[,2], mod5b[,2], mod7b[,2], mod9b[,2])

  ans <- matrix(ans, 1, length(ans))
  write.table(ans, file = "sim-t-mar-0101-result.txt", row.names = F, col.names = F, append = TRUE)
}

## write.table(result, file = "sim-t-mar-0101-result.txt", row.names = F, col.names = F)
sendEmail(subject="simulation-t-MAR", text="done", address="liuminzhao@gmail.com")

###############
## AFTER SIM; GET SUMMARY
###############
###############
## TRUE VALUE
###############

## quan1 <- function(y, x, tau){
##   return(tau - .5*pt((y - 2- x)/(1+alpha*x), df = 3) - .5*pt((y+2+x)/(1+alpha*x),df = 3))
## }

## SolveQuan1 <- function(x, tau){
##   uniroot(quan1, c(-30, 30), x = x, tau = tau)$root
## }

## xsim <- seq(0, 2, len = 100)
## y11 <- sapply(xsim, function(x) SolveQuan1(x, 0.1))
## y13 <- sapply(xsim, function(x) SolveQuan1(x, 0.3))
## y15 <- sapply(xsim, function(x) SolveQuan1(x, 0.5))
## y17 <- sapply(xsim, function(x) SolveQuan1(x, 0.7))
## y19 <- sapply(xsim, function(x) SolveQuan1(x, 0.9))

## q11 <- lm(y11~xsim)$coef
## q13 <- lm(y13~xsim)$coef
## q15 <- lm(y15~xsim)$coef
## q17 <- lm(y17~xsim)$coef
## q19 <- lm(y19~xsim)$coef

## quan2 <- function(y, x, tau){
##   return(tau - .5*pt((y  + 1.5*x)/(1+alpha*x), df = 3) - .5*pt((y - 2 + 0.5*x)/(1 + alpha*x),df = 3))
## }

## SolveQuan2 <- function(x, tau){
##   uniroot(quan2, c(-30, 30), x = x, tau = tau)$root
## }

## xsim <- seq(0, 2, len = 100)
## y25 <- sapply(xsim, function(x) SolveQuan2(x, 0.5))
## y29 <- sapply(xsim, function(x) SolveQuan2(x, 0.9))
## y27 <- sapply(xsim, function(x) SolveQuan2(x, 0.7))
## y23 <- sapply(xsim, function(x) SolveQuan2(x, 0.3))
## y21 <- sapply(xsim, function(x) SolveQuan2(x, 0.1))

## q21 <- lm(y21~xsim)$coef
## q23 <- lm(y23~xsim)$coef
## q25 <- lm(y25~xsim)$coef
## q27 <- lm(y27~xsim)$coef
## q29 <- lm(y29~xsim)$coef



## result <- read.table('sim-t-mar-0101-result.txt')
## trueq <- c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
## trueq <- rep(trueq, 3)

## ## MSE
## mse <- rep(0, 60)
## for (i in 1:60){
##   mse[i] <- mean((result[,i] - trueq[i])^2, trim = 0.05)
## }
## mseh2 <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))
## mserq <- rbind(matrix(mse[21:30], 2, 5), matrix(mse[31:40], 2, 5))
## msebb <- rbind(matrix(mse[41:50], 2, 5), matrix(mse[51:60], 2, 5))
## mytbl <- cbind(mseh2[,1], mserq[,1], msebb[,1], mseh2[,2], mserq[,2], msebb[,2], mseh2[,3], mserq[,3], msebb[,3], mseh2[,4], mserq[,4], msebb[,4], mseh2[,5], mserq[,5], msebb[,5])
## colnames(mytbl) <- rep(c('MM', 'RQ', 'BB'), 5)
## print(xtable(mytbl))
## print(mytbl)

## ## BIAS
## bias <- rep(0, 60)
## for (i in 1:60){
##   bias[i] <- mean((result[,i] - trueq[i]), trim = 0.05)
## }
## biash2 <- rbind(matrix(bias[1:10], 2, 5), matrix(bias[11:20], 2, 5))
## biasrq <- rbind(matrix(bias[21:30], 2, 5), matrix(bias[31:40], 2, 5))
## biasbb <- rbind(matrix(bias[41:50], 2, 5), matrix(bias[51:60], 2, 5))
## mytbl2 <- cbind(biash2[,1], biasrq[,1], biasbb[,1], biash2[,2], biasrq[,2], biasbb[,2], biash2[,3], biasrq[,3], biasbb[,3], biash2[,4], biasrq[,4], biasbb[,4], biash2[,5], biasrq[,5], biasbb[,5])
## colnames(mytbl2) <- rep(c('MM', 'RQ', 'BB'), 5)
## print(xtable(mytbl2))
## print(mytbl2)

## ## Efficiency
## efficiency <- apply(result[,1:60], 2, var)
## efficiencyh2 <- rbind(matrix(efficiency[1:10], 2, 5), matrix(efficiency[11:20], 2, 5))
## efficiencyrq <- rbind(matrix(efficiency[21:30], 2, 5), matrix(efficiency[31:40], 2, 5))
## efficiencybb <- rbind(matrix(efficiency[41:50], 2, 5), matrix(efficiency[51:60], 2, 5))
## mytbl3 <- cbind(efficiencyrq/efficiencyh2, efficiencybb/efficiencyh2)
## colnames(mytbl3) <- rep(c('RQ', 'BB'), 5)
## print(xtable(mytbl3))

## ## MCSE
## mcse <- rep(0, 60)
## for (i in 1:60){
##   mcse[i] <- sd((result[,i] - trueq[i])^2)/sqrt(boot)
## }
## mcseh2 <- rbind(matrix(mcse[1:10], 2, 5), matrix(mcse[11:20], 2, 5))
## mcserq <- rbind(matrix(mcse[21:30], 2, 5), matrix(mcse[31:40], 2, 5))
## mcsebb <- rbind(matrix(mcse[41:50], 2, 5), matrix(mcse[51:60], 2, 5))
## mytbl4 <- cbind(mcseh2[,1], mcserq[,1], mcsebb[,1], mcseh2[,2], mcserq[,2], mcsebb[,2], mcseh2[,3], mcserq[,3], mcsebb[,3], mcseh2[,4], mcserq[,4], mcsebb[,4], mcseh2[,5], mcserq[,5], mcsebb[,5])
## colnames(mytbl4) <- rep(c('MM', 'RQ', 'BB'), 5)
## print(xtable(mytbl4))
## print(mytbl4)

## ## combine mse and MCSE
## msemcse <- matrix(0, 4, 30)
## for (i in 1:15){
##   msemcse[, i*2 - 1] <- mytbl[, i]
##   msemcse[, i*2] <- mytbl4[, i]
## }
## print(msemcse)
## print(xtable(msemcse))

cat("Time: ", proc.time()[3] - start, '\n')

sink()
