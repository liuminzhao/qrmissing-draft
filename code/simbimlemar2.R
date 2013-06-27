## Time-stamp: <liuminzhao 04/17/2013 12:00:27>
## Simulation Bivariate case with MAR using heter2
rm(list = ls())
source('sendEmail.R')
source('BiMLESigma.R')
library(quantreg)
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
n <- 500
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

  mod1 <- BiQRGradient(y, R, x, tau = 0.1, method = "heter2")$param
  mod3 <- BiQRGradient(y, R, x, tau = 0.3, method = "heter2")$param
  mod5 <- BiQRGradient(y, R, x, tau = 0.5, method = "heter2")$param
  mod7 <- BiQRGradient(y, R, x, tau = 0.7, method = "heter2")$param
  mod9 <- BiQRGradient(y, R, x, tau = 0.9, method = "heter2")$param

  mod1rq <- as.vector(rq(y[,1]~x, tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)
  mod2rq <- as.vector(rq(y[,2][R==1]~x[R==1], tau = c(0.1, 0.3, 0.5, 0.7, 0.9))$coef)
  ans <- c(mod1[1:2], mod3[1:2], mod5[1:2], mod7[1:2], mod9[1:2],
           mod1[7:8], mod3[7:8], mod5[7:8], mod7[7:8], mod9[7:8], mod1rq, mod2rq)
##   ans <- c(mod1[1:2], mod3[1:2], mod5[1:2], mod7[1:2], mod9[1:2],
  ##         mod1[7:8], mod3[7:8], mod5[7:8], mod7[7:8], mod9[7:8])
  
}

save(result, file = "simbimar2Apr17.RData")
sendEmail(subject="simulation--bi-MAR2", text="done", address="liuminzhao@gmail.com")

trueq = c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
trueq <- rep(trueq, 2)
mse = rep(0, 40)
for (i in 1:40){
  mse[i] = mean((result[,i] - trueq[i])^2)
}
mse

## c(0.0676243087686589, 0.0419081503545594, 0.0442537968516964,
## 0.0496264730037774, 0.0270769491456292, 0.208023478367695, 0.0571549934633467,
## 0.0651301662591647, 0.0660934691145477, 0.0446580108037801, 0.0846528322932576,
## 0.0626178417798745, 0.0492628752831569, 0.0365766301832644, 0.0425691004993467,
## 0.0314464713910316, 0.0485128146259267, 0.0356912439098726, 0.0849878579549401,
## 0.0604820761972908)

## \begin{table}[ht]
## \begin{center}
## \begin{tabular}{rrrrrr}
##   \hline
##  & 1 & 2 & 3 & 4 & 5 \\ 
##   \hline
## 1 & 0.07 & 0.04 & 0.03 & 0.06 & 0.07 \\ 
##   2 & 0.04 & 0.05 & 0.21 & 0.07 & 0.04 \\ 
##   3 & 0.08 & 0.05 & 0.04 & 0.05 & 0.08 \\ 
##   4 & 0.06 & 0.04 & 0.03 & 0.04 & 0.06 \\ 
##    \hline
## \end{tabular}
## \end{center}
## \end{table}

## 04-17 Heter2 beta22
## > mseh2
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 0.07511 0.04440 0.02688 0.06220 0.05767
## [2,] 0.04834 0.03808 0.20632 0.07196 0.04132
## [3,] 0.09162 0.05654 0.04965 0.05561 0.07984
## [4,] 0.06760 0.04592 0.04038 0.04412 0.06434
## > mserq
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 0.09608 0.11119 0.26738 0.10873 0.09530
## [2,] 0.07151 0.08797 0.79691 0.09097 0.07543
## [3,] 0.11334 0.06537 0.05960 0.06496 0.11124
## [4,] 0.08002 0.05132 0.04685 0.05100 0.08596
