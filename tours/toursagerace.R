#!/bin/Rscript
## Time-stamp: <liuminzhao 06/06/2013 13:29:17>
## manipulate data TOURS
## 2013/06/05 focus on AGE and RACE
rm(list=ls())
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
weight2 <- weight2/100
weight3 <- weight3/100


n <- length(age)
y <- matrix(0, n, 2)
y[,1] <- weight2
y[,2] <- weight3


X <- matrix(0, n, 3)
X[,1] <- 1
X[,2] <- age_center
X[,3] <- race3

R <- 1 - as.numeric(is.na(weight3))

dat <- data.frame(weight2, weight3, trt, age_center, age=TOURS$AGE, race = factor(TOURS$RACE))

## PLOT
library(ggplot2)
library(gridExtra)
ggplot(data = dat, aes(x = age_center, y = weight2)) + geom_point()
ggplot(data = dat, aes(x = age_center, y = weight3)) + geom_point()

box1 <- ggplot(data = dat, aes(x = age, y = weight2*100)) + geom_point()+ ylab('Weight (Kg) after 6 months') + xlab('Age')
box2 <- ggplot(data = dat, aes(x = age, y = weight3*100)) + geom_point()+ ylab('Weight (Kg) after 18 months') + xlab('Age')
box3 <- ggplot(data = dat, aes(x = race, y = weight2*100)) + geom_boxplot() + scale_x_discrete(labels=c("Black", "White")) + ylab('Weight (Kg) after 6 months') + xlab('Race')
box4 <- ggplot(data = dat, aes(race, y = weight3*100)) + geom_boxplot() + ylab('Weight (Kg) after 18 months') + xlab('Race') + scale_x_discrete(labels=c("Black", "White"))
## scat1 <- ggplot(data = dat, aes(x = age_center, y = weight2, color = race)) + geom_point()
## scat2 <- ggplot(data = dat, aes(x = age_center, y = weight3, color = race)) + geom_point()
## scat3 <- ggplot(data = dat, aes(x = age_center, y = weight2, color = race)) + geom_point()
## scat4 <- ggplot(data = dat, aes(x = age_center, y = weight3, color = race)) + geom_point()

pdf('weight-plot-age.pdf')
sds <- grid.arrange(box1, box2, box3, box4, nrow = 2, ncol = 2)
dev.off()


###############
## ANALYSIS
###############

y[is.na(y[,2]),2] <- 0

source('~/Documents/qrmissing/code/BiMLESigma.R')
mod1 <- BiQRGradient(y, R, X, tau = 0.05, method = 'heter2')
mod3 <- BiQRGradient(y, R, X, tau = 0.3, method = 'heter2')
mod5 <- BiQRGradient(y, R, X, tau = 0.5, method = 'heter2')
mod7 <- BiQRGradient(y, R, X, tau = 0.7, method = 'heter2')
mod9 <- BiQRGradient(y, R, X, tau = 0.95, method = 'heter2')

coef1 <- mysummary(mod1)
coef3 <- mysummary(mod3)
coef5 <- mysummary(mod5)
coef7 <- mysummary(mod7)
coef9 <- mysummary(mod9)

coefw2 <- rbind(coef1[1,], coef3[1, ], coef5[1, ], coef7[1, ], coef9[1,])
coefw3 <- rbind(coef1[2,], coef3[2, ], coef5[2, ], coef7[2, ], coef9[2,])

rownames(coefw2) <- c('tau = 0.05', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.95')
rownames(coefw3) <- c('tau = 0.05', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.95')

colnames(coefw2) <- c('Intercept', 'Age(centered)', 'White')
colnames(coefw3) <- c('Intercept', 'Age(centered)', 'White')

library(xtable)
print(xtable(coefw2))
print(xtable(coefw3))

write.table(rbind(coefw2, coefw3), 'ageresult.txt', row.names=FALSE)

####################

## modrq = rq(change1 ~ age_center + factor(race), tau = seq(0.1,0.9,len =5))

## modrq2 = rq(change2[R==1] ~ age_center[R==1] + factor(race)[R==1], tau = seq(0.1,0.9,len =5))

## ## RESULTS

## result1 <- cbind(mod1$param[c(1:3, 13:15)], mod3$param[c(1:3, 13:15)],
##                  mod5$param[c(1:3, 13:15)], mod7$param[c(1:3, 13:15)],
##                  mod9$param[c(1:3, 13:15)])

## resultrq <- rbind(modrq$coef, modrq2$coef)

## library(xtable)
## print(xtable(result1))
## print(xtable(resultrq))
## write.table(result1, '../tours/result1.txt')
## write.table(resultrq, '../tours/resultrq.txt')
