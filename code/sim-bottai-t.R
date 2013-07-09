#!/bin/Rscript
##' Time-stamp: <liuminzhao 07/09/2013 13:09:33>
##' Simulation Bivariate case with MAR using heter2
##' Real MAR , not MCAR
##' correct heterogeneity parameters
##' test on Bottai's method
##' t distribution

sink('0709-bottai-t.txt')
rm(list = ls())
source('Bottai.R')
source('sendEmail.R')
library(quantreg)
library(xtable)
library(doMC)
registerDoMC()
options(cores = 8)
set.seed(3)

###############
## PARAMETER
###############
n <- 500
p <- 0.5
alpha <- 0.5

###############
## SIMULATION
###############
boot <- 1000
start <- proc.time()[3]

result <- foreach(icount(boot), .combine = rbind) %dopar% {
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

  mod1 <- Bottai(y, R, X, tau = 0.1)
  mod3 <- Bottai(y, R, X, tau = 0.3)
  mod5 <- Bottai(y, R, X, tau = 0.5)
  mod7 <- Bottai(y, R, X, tau = 0.7)
  mod9 <- Bottai(y, R, X, tau = 0.9)

  ans <- c(mod1[,1], mod3[,1], mod5[,1], mod7[,1], mod9[,1],
           mod1[,2], mod3[,2], mod5[,2], mod7[,2], mod9[,2])
}

write.table(result, file = "0709-bottai-t-result.txt", row.names = F, col.names = F)
sendEmail(subject="simulation-bi-bottai-t", text="done", address="liuminzhao@gmail.com")



###############
## TRUE VALUE
###############
quan1 <- function(y, x, tau){
  return(tau - .5*pt((y - 2- x)/(1+alpha*x), df = 3) - .5*pt((y+2+x)/(1+alpha*x),df = 3))
}

SolveQuan1 <- function(x, tau){
  uniroot(quan1, c(-30, 30), x = x, tau = tau)$root
}

xsim <- seq(0, 2, len = 100)
y11 <- sapply(xsim, function(x) SolveQuan1(x, 0.1))
y13 <- sapply(xsim, function(x) SolveQuan1(x, 0.3))
y15 <- sapply(xsim, function(x) SolveQuan1(x, 0.5))
y17 <- sapply(xsim, function(x) SolveQuan1(x, 0.7))
y19 <- sapply(xsim, function(x) SolveQuan1(x, 0.9))

q11 <- lm(y11~xsim)$coef
q13 <- lm(y13~xsim)$coef
q15 <- lm(y15~xsim)$coef
q17 <- lm(y17~xsim)$coef
q19 <- lm(y19~xsim)$coef

quan2 <- function(y, x, tau){
  return(tau - .5*pt((y  + 1.5*x)/(1+alpha*x), df = 3) - .5*pt((y - 2 + 0.5*x)/(1 + alpha*x),df = 3))
}

SolveQuan2 <- function(x, tau){
  uniroot(quan2, c(-30, 30), x = x, tau = tau)$root
}

xsim <- seq(0, 2, len = 100)
y25 <- sapply(xsim, function(x) SolveQuan2(x, 0.5))
y29 <- sapply(xsim, function(x) SolveQuan2(x, 0.9))
y27 <- sapply(xsim, function(x) SolveQuan2(x, 0.7))
y23 <- sapply(xsim, function(x) SolveQuan2(x, 0.3))
y21 <- sapply(xsim, function(x) SolveQuan2(x, 0.1))

q21 <- lm(y21~xsim)$coef
q23 <- lm(y23~xsim)$coef
q25 <- lm(y25~xsim)$coef
q27 <- lm(y27~xsim)$coef
q29 <- lm(y29~xsim)$coef

result <- read.table('0709-bottai-t-result.txt')
trueq <- c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
mse <- rep(0, 20)
for (i in 1:20){
  mse[i] <- mean((result[,i] - trueq[i])^2, trim = 0.05)
}

mse <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))

print(xtable(mse))

print(mse)
cat("Time: ", proc.time()[3] - start, '\n')
sink()
