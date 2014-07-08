## S2 Normal
rm(list = ls())
library(xtable)
alpha <- 0
quan1 <- function(y, x, tau){
  return(tau - .5*pnorm(y, 2+x, 1 + alpha*x) - .5*pnorm(y, -2-x, 1 + alpha*x))
}

SolveQuan1 <- function(x, tau){
  uniroot(quan1, c(-30, 30), x = x, tau = tau)$root
}

quan2 <- function(y, x, tau){
  return(tau - .5*pnorm(y, -1.5*x, (1+alpha*x)*sqrt(5/4)) - .5*pnorm(y, 4-.5*x, (1+alpha*x)*sqrt(5/4)))
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

result1 <- read.table('sim-normal-mnar-mar-0626-result.txt')

trueq <- c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
trueq <- rep(trueq, 4)
mse <- rep(0, 80)
for (i in 1:80){
  mse[i] <- mean((result1[,i] - trueq[i])^2)
}

msem1 <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))
msem2 <- rbind(matrix(mse[21:30], 2, 5), matrix(mse[31:40], 2, 5))
mserq <- rbind(matrix(mse[41:50], 2, 5), matrix(mse[51:60], 2, 5))
msebb <- rbind(matrix(mse[61:70], 2, 5), matrix(mse[71:80], 2, 5))

s2N <- cbind(c(msem1), c(msem2), c(mserq), c(msebb))
colnames(s2N) <- c('M1', 'M2', 'RQ', 'BB')
print(xtable(s2N))
print(s2N)

## s2 t3
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
  return(tau - .5*pt((y  + 1.5*x)/(1+alpha*x), df = 3) - .5*pt((y - 4 + 0.5*x)/(1 + alpha*x),df = 3))
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


result1 <- read.table('sim-t-mnar-mar-0626-result.txt')

trueq <- c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
trueq <- rep(trueq, 4)
mse <- rep(0, 80)
for (i in 1:80){
  mse[i] <- mean((result1[,i] - trueq[i])^2)
}

msem1 <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))
msem2 <- rbind(matrix(mse[21:30], 2, 5), matrix(mse[31:40], 2, 5))
mserq <- rbind(matrix(mse[41:50], 2, 5), matrix(mse[51:60], 2, 5))
msebb <- rbind(matrix(mse[61:70], 2, 5), matrix(mse[71:80], 2, 5))

s2t <- cbind(c(msem1), c(msem2), c(mserq), c(msebb))
colnames(s2t) <- c('M1', 'M2', 'RQ', 'BB')
print(xtable(s2t))
print(s2t)

## s2 LP
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
  return(quan - .5*pLD(y + (1.5*x), tau = tau/(1 + alpha*x)) - .5*pLD(y - (4 - 0.5*x),tau = tau/(1 + alpha*x)))
}

SolveQuan2 <- function(x, tau, quan){
  uniroot(quan2, c(-30, 30), x = x, tau = tau, quan = quan)$root
}

tau <- 1
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


result1 <- read.table('sim-lp-mnar-mar-0626-result.txt')

trueq <- c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
trueq <- rep(trueq, 4)
mse <- rep(0, 80)
for (i in 1:80){
  mse[i] <- mean((result1[,i] - trueq[i])^2)
}

msem1 <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))
msem2 <- rbind(matrix(mse[21:30], 2, 5), matrix(mse[31:40], 2, 5))
mserq <- rbind(matrix(mse[41:50], 2, 5), matrix(mse[51:60], 2, 5))
msebb <- rbind(matrix(mse[61:70], 2, 5), matrix(mse[71:80], 2, 5))

s2lp <- cbind(c(msem1), c(msem2), c(mserq), c(msebb))
colnames(s2lp) <- c('M1', 'M2', 'RQ', 'BB')
print(xtable(s2lp))
print(s2lp)

## s2 mix
quan1 <- function(y, x, tau){
  return(tau - .25*pnorm(y, 4+x, 1) - .25*pnorm(y, x, 1)
         - 0.25*pnorm(y, -x, 1) - 0.25*pnorm(y, -4-x, 1))
}

SolveQuan1 <- function(x, tau){
  uniroot(quan1, c(-30, 30), x = x, tau = tau)$root
}
quan2 <- function(y, x, tau){
  return(tau - 1/8*pnorm(y, -3 - 1.5*x, sqrt(5/4)) - 1/8*pnorm(y, 1 - 1.5*x, sqrt(5/4)) - 1/8*pnorm(y, -1 - 1.5*x, sqrt(5/4)) - 1/8*pnorm(y, 3 - 1.5*x, sqrt(5/4)) - 1/16*pnorm(y, -1 - .5*x, sqrt(5/4)) - 1/8*pnorm(y, 3 - .5*x, sqrt(5/4)) - 1/16*pnorm(y, 1 - .5*x, sqrt(5/4))- 1/8*pnorm(y, 5 - .5*x, sqrt(5/4)) - 1/16*pnorm(y, 5 - 1.5*x, sqrt(5/4)) - 1/16*pnorm(y, 7 - .5*x, sqrt(5/4)))
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

## result

result1 <- read.table('sim-m1-mnar-mar-result-0626.txt')

trueq <- c(q11, q13, q15, q17, q19, q21, q23, q25, q27, q29)
trueq <- rep(trueq, 4)
mse <- rep(0, 80)
for (i in 1:80){
  mse[i] <- mean((result1[,i] - trueq[i])^2)
}

msem1 <- rbind(matrix(mse[1:10], 2, 5), matrix(mse[11:20], 2, 5))
msem2 <- rbind(matrix(mse[21:30], 2, 5), matrix(mse[31:40], 2, 5))
mserq <- rbind(matrix(mse[41:50], 2, 5), matrix(mse[51:60], 2, 5))
msebb <- rbind(matrix(mse[61:70], 2, 5), matrix(mse[71:80], 2, 5))

s2mix <- cbind(c(msem1), c(msem2), c(mserq), c(msebb))
colnames(s2mix) <- c('M1', 'M2', 'RQ', 'BB')
print(xtable(s2mix))
print(s2mix)

## merge

s2tbl <- cbind(s2N, s2t, s2lp, s2mix)
print(xtable(s2tbl))
