#!/bin/Rscript
##' Time-stamp: <liuminzhao 07/07/2013 14:54:10>
##' Simulation Bivariate case with MAR using heter2
##' Real MAR , not MCAR
##' correct heterogeneity parameters

sink('0705.txt')
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
## PARAMETER
###############
n <- 500
p <- 0.5
alpha <- 0.5

###############
## SIMULATION
###############
boot <- 8
count <- rep(0, 5)
start <- proc.time()[3]

result <- foreach(icount(boot), .combine = rbind) %dopar% {
  R <- rbinom(n, 1, p)
  x <- runif(n, 0, 2)
  y <- matrix(0, n, 2)
  for (i in 1:n){
    if (R[i] == 1){
      y[i, 1] <- 2 + x[i] +(1 + 0.5*x[i])*rnorm(1)
      y[i, 2] <- 1 - x[i] - 0.5 * y[i, 1] + (1+0.5*x[i])*rnorm(1)
    } else {
      y[i, 1] <- -2 - x[i] +(1 + 0.5*x[i])*rnorm(1)
      y[i, 2] <- 1 - x[i] - 0.5 * y[i, 1] + (1+0.5*x[i])*rnorm(1)
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

  count[1] <- count[1] + mod1$converge
  count[2] <- count[2] + mod3$converge
  count[3] <- count[3] + mod5$converge
  count[4] <- count[4] + mod7$converge
  count[5] <- count[5] + mod9$converge

  mod1rq <- as.vector(rq(y[,1]~x, tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)
  mod2rq <- as.vector(rq(y[,2][R==1]~x[R==1], tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)

  ans <- c(mod1mm[1,], mod3mm[1,], mod5mm[1,], mod7mm[1,], mod9mm[1,],
           mod1mm[2,], mod3mm[2,], mod5mm[2,], mod7mm[2,], mod9mm[2,], mod1rq, mod2rq)
}

write.table(result, file = "0705-result.txt", row.names = F, col.names = F)
sendEmail(subject="simulation-bi--MAR", text="done", address="liuminzhao@gmail.com")

###############
## TRUE VALUE
###############
quan1 <- function(y, x, tau){
  return(tau - .5*pnorm(y, 2+x, 1 + 0.5*x) - .5*pnorm(y, -2-x, 1 + 0.5*x))
}

SolveQuan1 <- function(x, tau){
  uniroot(quan1, c(-30, 30), x = x, tau = tau)$root
}

quan2 <- function(y, x, tau){
  return(tau - .5*pnorm(y, -1.5*x, (1+0.5*x)*sqrt(5/4)) - .5*pnorm(y, 2-.5*x, (1+0.5*x)*sqrt(5/4)))
}

SolveQuan2 <- function(x, tau){
  uniroot(quan2, c(-30, 30), x = x, tau = tau)$root
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

result <- read.table('0705-result.txt')
trueq <- c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
trueq <- rep(trueq, 2)
mse <- rep(0, 40)
for (i in 1:40){
  mse[i] <- mean((result[,i] - trueq[i])^2, trim = 0.05)
}

mseh2 <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))
mserq <- rbind(matrix(mse[21:30], 2, 5), matrix(mse[31:40], 2, 5))
mytbl <- cbind(mseh2[,1], mserq[,1], mseh2[,2], mserq[,2], mseh2[,3], mserq[,3], mseh2[,4], mserq[,4], mseh2[,5], mserq[,5])

print(xtable(mseh2))
print(xtable(mserq))

colnames(mytbl) <- rep(c('MM', 'RQ'), 5)

print(xtable(mytbl))

cat("Time: ", proc.time()[3] - start, '\n')
cat("Converged: ", count)
sink()
