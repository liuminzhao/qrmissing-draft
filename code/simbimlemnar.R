## Time-stamp: <liuminzhao 03/29/2013 21:18:45>
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

quan2 <- function(y, x, tau){
  return(tau - .5*pnorm(y, 1-x, 1) - .5*pnorm(y, 2-x,1))
}

SolveQuan2 <- function(x, tau){
  uniroot(quan2, c(-30, 30), x = x, tau = tau)$root
}

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

  mod1 <- BiQRGradient(y, R, x, tau = 0.1, sp = c(1,0,0,1))$param
  mod3 <- BiQRGradient(y, R, x, tau = 0.3, sp = c(1,0,0,1))$param
  mod5 <- BiQRGradient(y, R, x, tau = 0.5, sp = c(1,0,0,1))$param
  mod7 <- BiQRGradient(y, R, x, tau = 0.7, sp = c(1,0,0,1))$param
  mod9 <- BiQRGradient(y, R, x, tau = 0.9, sp = c(1,0,0,1))$param

  ans <- c(mod1[1:2], mod3[1:2], mod5[1:2], mod7[1:2], mod9[1:2],
           mod1[7:8], mod3[7:8], mod5[7:8], mod7[7:8], mod9[7:8])
  
}

save(result, file = "simbimnar.RData")
sendEmail(subject="simulation--b1-MNAR", text="done", address="liuminzhao@gmail.com")

trueq = c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
mse = rep(0, 20)
for (i in 1:20){
  mse[i] = mean((result[,i] - trueq[i])^2)
}
mse


