## Time-stamp: <liuminzhao 03/29/2013 15:54:31>
## Simulation Bivariate case with MAR
rm(list = ls())
source('sendEmail.R')
source('BiMLESigma.R')
library(doMC)
registerDoMC()
options(cores = 8)
set.seed(1)


###############
## TRUE VALUE 
###############

quan1 <- function(y, x, tau){
  return(tau - .5*pnorm(y, 1+x, 1) - .5*pnorm(y, -1-x,1))
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

q21 <- c(1 +  qnorm(0.1), -1)
q23 <- c(1 +  qnorm(0.3), -1)
q25 <- c(1 +  qnorm(0.5), -1)
q27 <- c(1 +  qnorm(0.7), -1)
q29 <- c(1 +  qnorm(0.9), -1)


###############
## PARAMETER 
###############
p <- 0.5
n <- 200
b01 <- 1
b00 <- -1
b11 <- 1
b10 <- -1
sigma1 <- 1
sigma0 <- 1
b02 <- 1
b12 <- -1
sigma2 <- 1


###############
## SIMULATION 
###############
boot <- 1000
resultsave <- matrix(0, boot, 20)

result <- foreach(icount(boot), .combine = rbind) %dopar% {
  R <- rbinom(n, 1, p)
  x <- runif(n, 0, 2)
  y <- matrix(0, n, 2)
  for (i in 1:n){
    if (R[i] == 1){
      y[i, 1] <- rnorm(1, b01 + x[i] * b11, sigma1)
      y[i, 2] <- rnorm(1, b02 + x[i] * b12, sigma2)
    } else {
      y[i, 1] <- rnorm(1, b00 + x[i] * b10, sigma0)
      y[i, 2] <- rnorm(1, b00 + x[i] * b10, sigma0)
    }
  }

  mod1 <- BiQRGradient(y, R, x, tau = 0.1)$param
  mod3 <- BiQRGradient(y, R, x, tau = 0.3)$param
  mod5 <- BiQRGradient(y, R, x, tau = 0.5)$param
  mod7 <- BiQRGradient(y, R, x, tau = 0.7)$param
  mod9 <- BiQRGradient(y, R, x, tau = 0.9)$param

  ans <- c(mod1[1:2], mod3[1:2], mod5[1:2], mod7[1:2], mod9[1:2],
           mod1[7:8], mod3[7:8], mod5[7:8], mod7[7:8], mod9[7:8])
  
}

save(result, file = "simbimar.RData")
sendEmail(subject="simulation--b1-MAR", text="done", address="liuminzhao@gmail.com")
