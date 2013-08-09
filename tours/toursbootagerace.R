#!/bin/Rscript
## Time-stamp: <liuminzhao 08/08/2013 23:02:00>
## bootstrap on tours data
## weight2 and weight3
## scaled by 1/100
## covariates: 3 treatments and 2 races
## 2013/07/29 add init1 and init9 , modfiied saved file name, change to 10 core
## 2013/08/03 Using QRMissingBi.R and change boot to 50
## 2013/08/08 on Int, age, Race
rm(list=ls())
library(compiler)
enableJIT(3)
enableJIT(3)
library(rootSolve)
library(quantreg)
library(minqa)
source('../code/QRMissingBi.R')
source('~/Documents/qrmissing/code/sendEmail.R')
library(doMC)
registerDoMC()
options(cores = 10)
set.seed(1)

TOURS <- read.csv('~/Documents/qrmissing/tours/tours.csv')
TOURS <- subset(TOURS, RACE==1 | RACE==3)
weight1 <- TOURS$wtkg1
weight2 <- TOURS$wtkg2
weight3 <- TOURS$wtkg3
age <- TOURS$AGE
trt <- TOURS$TREATMENT
race3 <- as.numeric(TOURS$RACE == 3)
age_center <- (age-mean(age))/sd(age)

## scaled by 1/100
weight1 <- weight1/100
weight2 <- weight2/100
weight3 <- weight3/100

## new dataset
NEWTOURS <- data.frame(w2 = weight2, w3 = weight3, age_center, race3, w1 = weight1)

###############
## BOOTSTRAP
###############
boot <- 10
n <- dim(NEWTOURS)[1]
y <- matrix(0, n, 2)
X <- matrix(0, n, 3)
X[,1] <- 1

start <- proc.time()[3]

QRMissingBic <- cmpfun(QRMissingBi)

result <- foreach(icount(boot), .combine = rbind) %dopar% {
  indices <- sample(n, replace = TRUE)
  dat <- NEWTOURS[indices, ]
  y[, 1] <- dat$w2
  y[, 2] <- dat$w3
  X[,2] <- dat$age_center
  X[,3] <- dat$race3
  X[,4] <- dat$w1
  R <- 1 - as.numeric(is.na(dat$w3))
  y[is.na(y[,2]),2] <- 0

  mod1 <- QRMissingBic(y, R, X, tau = 0.1)
  mod3 <- QRMissingBic(y, R, X, tau = 0.3)
  mod5 <- QRMissingBic(y, R, X, tau = 0.5)
  mod7 <- QRMissingBic(y, R, X, tau = 0.7)
  mod9 <- QRMissingBic(y, R, X, tau = 0.9)

  coef1 <- coef(mod1)
  coef3 <- coef(mod3)
  coef5 <- coef(mod5)
  coef7 <- coef(mod7)
  coef9 <- coef(mod9)
  count1 <- mod1$ierr
  count3 <- mod3$ierr
  count5 <- mod5$ierr
  count7 <- mod7$ierr
  count9 <- mod9$ierr

  coefw2 <- c(coef1[1,], coef3[1, ], coef5[1, ], coef7[1, ], coef9[1,])
  coefw3 <- c(coef1[2,], coef3[2, ], coef5[2, ], coef7[2, ], coef9[2,])

  ans <- c(coefw2, coefw3, count1, count3, count5, count7, count9)

}

write.table(result, file = "toursbootagerace-0808.txt", row.names = FALSE, col.names = FALSE)
sendEmail(subject="boot tours base", text="done", address="liuminzhao@gmail.com")

print(proc.time()[3] - start)
