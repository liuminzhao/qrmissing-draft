## Time-stamp: <liuminzhao 04/17/2013 12:51:08>
## Simulation Bivariate case with MAR using heter2
## Real MAR , not MCAR

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
quan2 <- function(y, x, tau){
  return(tau - .5*pnorm(y, 2+x, 1 + 0.5*x) - .5*pnorm(y, -2-x, 1 + 0.5*x))
}

SolveQuan2 <- function(x, tau){
  uniroot(quan2, c(-30, 30), x = x, tau = tau)$root
}

quan3 <- function(y, x, tau){
  return(tau - .5*pnorm(y, -1.5*x, 5/4 + 1/8*x) - .5*pnorm(y, 2-.5*x, 5/4 + 1/8*x))
}

SolveQuan3 <- function(x, tau){
  uniroot(quan3, c(-30, 30), x = x, tau = tau)$root
}

xx <- seq(0, 2, len = 100)
quan <- seq(0.1, 0.9, len= 5)
qq3 <- sapply(quan, function(t) sapply(xx, function(x) SolveQuan3(x, t)))
(qest3 <- sapply(quan, function(t) lm(sapply(xx, function(x) SolveQuan3(x, t)) ~ xx)$coef))
qq2 <- sapply(quan, function(t) sapply(xx, function(x) SolveQuan2(x, t)))
(qest2 <- sapply(quan, function(t) lm(sapply(xx, function(x) SolveQuan2(x, t)) ~ xx)$coef))



###############
## PARAMETER 
###############
n <- 500
p <- 0.5


###############
## SIMULATION 
###############
boot <- 100

result <- foreach(icount(boot), .combine = rbind) %dopar% {
  R <- rbinom(n, 1, p)
  x <- runif(n, 0, 2)
  y <- matrix(0, n, 2)
  for (i in 1:n){
    if (R[i] == 1){
      y[i, 1] <- 2 + x[i] +(1 + 0.5*x[i])*rnorm(1)
      y[i, 2] <- 1 - x[i] - 0.5 * y[i, 1] + rnorm(1)
    } else {
      y[i, 1] <- -2 - x[i] +(1 + 0.5*x[i])*rnorm(1)
      y[i, 2] <- 1 - x[i] - 0.5 * y[i, 1] + rnorm(1)
    }
  }

  mod1 <- BiQRGradient(y, R, x, tau = 0.1, method = 'heter2')$param
  mod3 <- BiQRGradient(y, R, x, tau = 0.3, method = 'heter2')$param
  mod5 <- BiQRGradient(y, R, x, tau = 0.5, method = 'heter2')$param
  mod7 <- BiQRGradient(y, R, x, tau = 0.7, method = 'heter2')$param
  mod9 <- BiQRGradient(y, R, x, tau = 0.9, method = 'heter2')$param

  mod1rq <- as.vector(rq(y[,1]~x, tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)
  mod2rq <- as.vector(rq(y[,2][R==1]~x[R==1], tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)
  ans <- c(mod1[1:2], mod3[1:2], mod5[1:2], mod7[1:2], mod9[1:2],
           mod1[7:8], mod3[7:8], mod5[7:8], mod7[7:8], mod9[7:8], mod1rq, mod2rq)
}

save(result, file = "0417h2.RData")
sendEmail(subject="simulation--bi-Heter2-MAR", text="done", address="liuminzhao@gmail.com")

trueq <- c(c(qest2), c(qest3))
trueq <- rep(trueq, 2)
mse = rep(0, 40)
for (i in 1:40){
  mse[i] = mean((result[,i] - trueq[i])^2)
}

mseh2 <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))
mserq <- rbind(matrix(mse[21:30], 2, 5), matrix(mse[31:40], 2, 5))
print(mseh2)
print(mserq)
print(xtable(mseh2))
print(xtable(mserq))

## R > print(mseh2)
##            [,1]       [,2]       [,3]       [,4]       [,5]
## [1,] 0.08603033 0.11830813 0.10919735 0.16449038 0.10429721
## [2,] 0.08955918 0.06691685 0.13636475 0.07675017 0.09921015
## [3,] 0.07656087 0.06661742 0.06497826 0.12139400 0.24360734
## [4,] 0.06190360 0.04933893 0.05722029 0.06961475 0.08707117
## R > print(mserq)
##           [,1]      [,2]      [,3]      [,4]      [,5]
## [1,] 0.1475245 0.1931717 1.0750455 0.1880457 0.1452128
## [2,] 0.1512498 0.1921511 1.1871157 0.1983399 0.1519454
## [3,] 0.2737565 0.5946287 1.0826260 1.7478138 2.9221457
## [4,] 0.1680383 0.1253706 0.3297871 0.7473080 0.9648093
