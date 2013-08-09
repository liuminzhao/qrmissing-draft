#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/09/2013 14:05:59>
##' manipulate data TOURS
##' 2013/06/05 focus on AGE and RACE
##' 2013/06/22 add baseline y0 as a covariate
##' 2013/08/02 using QRMissingBi

rm(list=ls())
library(compiler)
enableJIT(3)
enableJIT(3)
library(rootSolve)
library(quantreg)
library(minqa)
library(numDeriv)
source('../code/QRMissingBi.R')

TOURS <- read.csv('~/Documents/qrmissing/tours/tours.csv')

TOURS <- subset(TOURS, RACE==1 | RACE==3)

weight1 <- TOURS$wtkg1
weight2 <- TOURS$wtkg2
weight3 <- TOURS$wtkg3
age <- TOURS$AGE
trt <- TOURS$TREATMENT
age_center <- (age-mean(age))/sd(age)
race3 <- as.numeric(TOURS$RACE == 3)

## center weight2?
## weight2 <- scale(weight2)
## weight3 <- scale(weight3)
weight1 <- weight1/100
weight2 <- weight2/100
weight3 <- weight3/100


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

## PLOT
## library(ggplot2)
## library(gridExtra)
## ggplot(data = dat, aes(x = age_center, y = weight2)) + geom_point()
## ggplot(data = dat, aes(x = age_center, y = weight3)) + geom_point()

## box1 <- ggplot(data = dat, aes(x = age, y = weight2*100)) + geom_point()+ ylab('Weight (Kg) at 6 months') + xlab('Age')
## box2 <- ggplot(data = dat, aes(x = age, y = weight3*100)) + geom_point()+ ylab('Weight (Kg) at 18 months') + xlab('Age')
## box3 <- ggplot(data = dat, aes(x = race, y = weight2*100)) + geom_boxplot() + scale_x_discrete(labels=c("Black", "White")) + ylab('Weight (Kg) at 6 months') + xlab('Race')
## box4 <- ggplot(data = dat, aes(race, y = weight3*100)) + geom_boxplot() + ylab('Weight (Kg) at 18 months') + xlab('Race') + scale_x_discrete(labels=c("Black", "White"))
## box5 <- ggplot(data = dat, aes(x = weight1*100, y = weight2*100)) + geom_point() + ylab('Weight (Kg) at 6 months') + xlab('Weight (Kg) at baseline')
## box6 <- ggplot(data = dat, aes(x = weight1*100, y = weight3*100)) + geom_point() + ylab('Weight (Kg) at 18 months') + xlab('Weight (Kg) at baseline')

## scat1 <- ggplot(data = dat, aes(x = age_center, y = weight2, color = race)) + geom_point()
## scat2 <- ggplot(data = dat, aes(x = age_center, y = weight3, color = race)) + geom_point()
## scat3 <- ggplot(data = dat, aes(x = age_center, y = weight2, color = race)) + geom_point()
## scat4 <- ggplot(data = dat, aes(x = age_center, y = weight3, color = race)) + geom_point()

## pdf('weight-age-race-base.pdf')
## sds <- grid.arrange(box1, box2, box3, box4, box5, box6, nrow = 2, ncol = 3)
## dev.off()


###############
## ANALYSIS
###############

## y[is.na(y[,2]),2] <- 0
QRMissingBic <- cmpfun(QRMissingBi)

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

coefw2 <- rbind(coef1[1,], coef3[1, ], coef5[1, ], coef7[1, ], coef9[1,])
coefw3 <- rbind(coef1[2,], coef3[2, ], coef5[2, ], coef7[2, ], coef9[2,])

rownames(coefw2) <- c('tau = 0.1', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.9')
rownames(coefw3) <- c('tau = 0.1', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.9')

colnames(coefw2) <- c('Intercept', 'Age(centered)', 'White', 'BaseWeight')
colnames(coefw3) <- c('Intercept', 'Age(centered)', 'White', 'BaseWeight')

## sd
sd1 <- mod1$se
sd3 <- mod3$se
sd5 <- mod5$se
sd7 <- mod7$se
sd9 <- mod9$se

sdw2 <- rbind(sd1[1,], sd3[1, ], sd5[1, ], sd7[1, ], sd9[1,])
sdw3 <- rbind(sd1[2,], sd3[2, ], sd5[2, ], sd7[2, ], sd9[2,])

rownames(sdw2) <- c('tau = 0.1', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.9')
rownames(sdw3) <- c('tau = 0.1', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.9')

colnames(sdw2) <- c('Intercept', 'Age(centered)', 'White')
colnames(sdw3) <- c('Intercept', 'Age(centered)', 'White')

## merge
w2 <- cbind(coefw2[,1], sdw2[,1],coefw2[,2], sdw2[,2], coefw2[,3], sdw2[,3], coefw2[,4], sdw2[,4])
w3 <- cbind(coefw3[,1], sdw3[,1],coefw3[,2], sdw3[,2], coefw3[,3], sdw3[,3], coefw3[,4], sdw3[,4])
colnames(w2) <- c('Intercept', 'sd', 'Age(centered)', 'sd', 'White', 'sd', 'Baseweight', 'sd')
colnames(w3) <- c('Intercept', 'sd', 'Age(centered)', 'sd', 'White', 'sd', 'Baseweight', 'sd')

library(xtable)
print(xtable(w2))
print(xtable(w3))

write.table(rbind(w2, w3), 'ageracebase.txt', row.names=FALSE)
