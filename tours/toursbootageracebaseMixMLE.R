#!/bin/Rscript
## Time-stamp: <liuminzhao 05/10/2014 12:02:01>
## bootstrap on tours data
## weight2 and weight3
## scaled by 1/100
## covariates: 3 treatments and 2 races
## 2013/07/29 add init1 and init9 , modfiied saved file name, change to 10 core
## 2013/08/03 Using QRMissingBi.R and change boot to 50
## 2013/08/08 on Int, age, Race , base
##' 2013/08/26 using qrmissing package
rm(list=ls())
library(qrmissing)
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
age_center <- (age-50)/5

## scaled by 1/100
weight1 <- weight1/10
weight2 <- weight2/10
weight3 <- weight3/10

## new dataset
NEWTOURS <- data.frame(w2 = weight2, w3 = weight3, age_center, race3, w1 = weight1)

###############
## BOOTSTRAP
###############
boot <- 1000
n <- dim(NEWTOURS)[1]
y <- matrix(0, n, 2)
X <- matrix(0, n, 4)
X[,1] <- 1

start <- proc.time()[3]

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

    mod1 <- QRMissingBiMixMLE(y, R, X, tau = 0.1, K = 2)
    mod3 <- QRMissingBiMixMLE(y, R, X, tau = 0.3, K = 2)
    mod5 <- QRMissingBiMixMLE(y, R, X, tau = 0.5, K = 2)
    mod7 <- QRMissingBiMixMLE(y, R, X, tau = 0.7, K = 2)
    mod9 <- QRMissingBiMixMLE(y, R, X, tau = 0.9, K = 2)

    coef1 <- rbind(coef(mod1)$gamma1, coef(mod1)$gamma2)
    coef3 <- rbind(coef(mod3)$gamma1, coef(mod3)$gamma2)
    coef5 <- rbind(coef(mod5)$gamma1, coef(mod5)$gamma2)
    coef7 <- rbind(coef(mod7)$gamma1, coef(mod7)$gamma2)
    coef9 <- rbind(coef(mod9)$gamma1, coef(mod9)$gamma2)

    count1 <- mod1$ierr
    count3 <- mod3$ierr
    count5 <- mod5$ierr
    count7 <- mod7$ierr
    count9 <- mod9$ierr

    coef1[, c(1, 2, 3)] <- 10 * coef1[, c(1, 2, 3)]
    coef3[, c(1, 2, 3)] <- 10 * coef3[, c(1, 2, 3)]
    coef5[, c(1, 2, 3)] <- 10 * coef5[, c(1, 2, 3)]
    coef7[, c(1, 2, 3)] <- 10 * coef7[, c(1, 2, 3)]
    coef9[, c(1, 2, 3)] <- 10 * coef9[, c(1, 2, 3)]

    coefw2 <- c(coef1[1,], coef3[1, ], coef5[1, ], coef7[1, ], coef9[1,])
    coefw3 <- c(coef1[2,], coef3[2, ], coef5[2, ], coef7[2, ], coef9[2,])

    ans <- c(coefw2, coefw3, count1, count3, count5, count7, count9)

}

write.table(result, file = "toursbootageracebaseMixMLE-0510.txt", row.names = FALSE, col.names = FALSE)

print(proc.time()[3] - start)
