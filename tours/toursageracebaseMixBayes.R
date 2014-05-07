#!/bin/Rscript
##' Time-stamp: <liuminzhao 05/07/2014 16:29:59>
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

mcmc <- list(nsave = 10000, nskip = 2, nburn = 4000, ndisp = 1000)
## prior
q <- dim(X)[2]
betapm <- 0
gammapm <- rep(0, q)
betapv <- 1
gammapv <- rep(10, q)
beta2pm <- 0
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

mod1 <- QRMissingBiBayesMix(y, R, X, tau = 0.1, mcmc = mcmc, prior = prior, method = "DP", sampling = 'element')
mod3 <- QRMissingBiBayesMix(y, R, X, tau = 0.3, mcmc = mcmc, prior = prior, method = "DP", sampling = 'element')
mod5 <- QRMissingBiBayesMix(y, R, X, tau = 0.5, mcmc = mcmc, prior = prior, method = "DP", sampling = 'element')
mod7 <- QRMissingBiBayesMix(y, R, X, tau = 0.7, mcmc = mcmc, prior = prior, method = "DP", sampling = 'element')
mod9 <- QRMissingBiBayesMix(y, R, X, tau = 0.9, mcmc = mcmc, prior = prior, method = "DP", sampling = 'element')

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

write.table(rbind(coefw2, coefw3), 'ageracebaseMixMLE.txt', row.names=FALSE)
