#!/bin/Rscript
##' Time-stamp: <liuminzhao 05/11/2014 11:15:59>
##' manipulate data TOURS
##' 2013/06/05 focus on AGE and RACE
##' 2013/06/22 add baseline y0 as a covariate
##' 2013/08/02 using QRMissingBi
##' 2013/08/25 using qrmissing package
##' 2014-05-07 using Bayes Mix K = 2

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
mcmc <- list(nsave = 30000, nskip = 2, nburn = 10000, ndisp = 10000)

## prior
q <- dim(X)[2]
betapm <- 0
gammapm <- rep(0, q)
betapv <- 1
gammapv <- rep(10, q)
beta2pm <- 0.36
beta2pv <- 0.001
K <- 2
betaysppm <- sigma21sppm <- 0
betaysppv <- sigma21sppv <- 1

prior <- list(betapm = betapm, betapv = betapv,
               gammapm = gammapm,
               gammapv= gammapv,
               beta2pm = beta2pm, beta2pv = beta2pv,
               alpha1 = 1, alpha2 = 1,
               K = K,
               sigmamu = 2,
               sigmaa = 2, sigmab = 2,
               eta1 = 2, eta2 = 2,
               A = sqrt(1000),
               tunegamma1 = c(0.07, 0.05, 0.05, 0.07),
               tunegamma2 = c(0.07, 0.05, 0.05, 0.07),
               tunebeta2sp = 0.001,
               tunebeta1 = 0.03,
               tunebetay = 0.03,
               tunesigma = 0.03,
               tunep = 0.03,
               arate = 0.25,
               tunemu = 0.03,
               tuneomega = 30)

## y[is.na(y[,2]),2] <- 0

mod9 <- ToursMNARBayesMix(y, R, X, tau = 0.9, mcmc = mcmc, prior = prior, method = "DP", sampling = 'element')

coef9 <- rbind(coef(mod9)$gamma1, coef(mod9)$gamma2)

coef9[, c(1, 2, 3)] <- 10 * coef9[, c(1, 2, 3)]

coefw2 <- coef9[1,]
coefw3 <- coef9[2,]

arrange <- function(mod){
    cigamma1 <- confintQRMissingBiBayesMix(mod)$gamma1
    cigamma1[, 1:3] <- 10 * cigamma1[, 1:3]
    cigamma1 <- c(cigamma1)

    cigamma2 <- confintQRMissingBiBayesMix(mod)$gamma2
    cigamma2[, 1:3] <- 10 * cigamma2[, 1:3]
    cigamma2 <- c(cigamma2)

    return(rbind(cigamma1, cigamma2))

}

ci9 <- arrange(mod9)

ciw2 <- ci9[1, ]
ciw3 <- ci9[2, ]

matw2 <- c(coefw2, ciw2)
matw2 <- matw2[c(1, 5,6, 2, 7, 8, 3, 9, 10, 4, 11, 12)]

matw3 <- c(coefw3, ciw3)
matw3 <- matw3[c(1, 5,6, 2, 7, 8, 3, 9, 10, 4, 11, 12)]

write.table(rbind(matw2, matw3), 'ageracebaseMixBayes-mnar-0511-9.txt', row.names=FALSE)
