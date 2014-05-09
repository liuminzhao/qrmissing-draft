#!/bin/Rscript
##' Time-stamp: <liuminzhao 05/09/2014 10:10:05>
##' 2013/06/05 focus on AGE and RACE
##' 2013/06/22 add baseline y0 as a covariate
##' 2013/07/05 MNAR
##' 2013/08/02 Using ToursMNAR2.R (R optimization of QRMissingBi.R)
##' 2013/08/07 Using ToursMNAR3.R (uobyqa optimization and 0.036 + Y1 for Y2|Y1, R = 0
##' 2013/08/26 using qrmissing package
##' 2014-05-07 use ToursMNARMixMLE
rm(list=ls())
library(qrmissing)

TOURS <- read.csv('~/Documents/qrmissing/tours/tours.csv')

TOURS <- subset(TOURS, RACE==1 | RACE==3)

weight1 <- TOURS$wtkg1
weight2 <- TOURS$wtkg2
weight3 <- TOURS$wtkg3
age <- TOURS$AGE
trt <- TOURS$TREATMENT
age_center <- (age-50)/5
race3 <- as.numeric(TOURS$RACE == 3)

## center weight2?
## weight2 <- scale(weight2)
## weight3 <- scale(weight3)
weight1 <- weight1/10
weight2 <- weight2/10
weight3 <- weight3/10

n <- length(age)
y <- matrix(0, n, 2)
y[,1] <- weight2
y[,2] <- weight3


X <- matrix(0, n, 4)
X[,1] <- 1
X[,2] <- age_center
X[,3] <- race3
X[,4] <- weight1

R <- 1 - as.numeric(is.na(weight3))

dat <- data.frame(weight2, weight3, trt, age_center, age=TOURS$AGE, race = factor(TOURS$RACE), base = weight1)

###############
## ANALYSIS
###############

## y[is.na(y[,2]),2] <- 0


mod1 <- ToursMNARMixMLE(y, R, X, tau = 0.1, K = 2)
mod3 <- ToursMNARMixMLE(y, R, X, tau = 0.3, K = 2)
mod5 <- ToursMNARMixMLE(y, R, X, tau = 0.5, K = 2)
mod7 <- ToursMNARMixMLE(y, R, X, tau = 0.7, K = 2)
mod9 <- ToursMNARMixMLE(y, R, X, tau = 0.9, K = 2)

print(mod1$ierr)
print(mod3$ierr)
print(mod5$ierr)
print(mod7$ierr)
print(mod9$ierr)

  coef1 <- rbind(coef(mod1)$gamma1, coef(mod1)$gamma2)
  coef3 <- rbind(coef(mod3)$gamma1, coef(mod3)$gamma2)
  coef5 <- rbind(coef(mod5)$gamma1, coef(mod5)$gamma2)
  coef7 <- rbind(coef(mod7)$gamma1, coef(mod7)$gamma2)
  coef9 <- rbind(coef(mod9)$gamma1, coef(mod9)$gamma2)


coef1[, c(1, 2, 3)] <- 10 * coef1[, c(1, 2, 3)]
coef3[, c(1, 2, 3)] <- 10 * coef3[, c(1, 2, 3)]
coef5[, c(1, 2, 3)] <- 10 * coef5[, c(1, 2, 3)]
coef7[, c(1, 2, 3)] <- 10 * coef7[, c(1, 2, 3)]
coef9[, c(1, 2, 3)] <- 10 * coef9[, c(1, 2, 3)]

coefw2 <- rbind(coef1[1,], coef3[1, ], coef5[1, ], coef7[1, ], coef9[1,])
coefw3 <- rbind(coef1[2,], coef3[2, ], coef5[2, ], coef7[2, ], coef9[2,])

rownames(coefw2) <- c('tau = 0.1', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.9')
rownames(coefw3) <- c('tau = 0.1', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.9')

colnames(coefw2) <- c('Intercept', 'Age(centered)', 'White', 'BaseWeight')
colnames(coefw3) <- c('Intercept', 'Age(centered)', 'White', 'BaseWeight')

library(xtable)
print(xtable(coefw2))
print(xtable(coefw3))

write.table(rbind(coefw2, coefw3), 'ageracebasemnar-mix-mle-0508.txt', row.names=FALSE)
