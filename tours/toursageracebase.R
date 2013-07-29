#!/bin/Rscript
##' Time-stamp: <liuminzhao 07/29/2013 04:04:00>
##' manipulate data TOURS
##' 2013/06/05 focus on AGE and RACE
##' 2013/06/22 add baseline y0 as a covariate

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

y[is.na(y[,2]),2] <- 0

## initial values
init9 <- c(0.032587903868868, 0.00490885508829561, -0.0260610021043252,
0.965537273955055, 0.0422616997543943, 0.00332924248002095, -0.0228402335710618,
-0.0247757293042863, -3.89490034005928, 0.00105545973672419,
0.240487850131395, 0.672555884353915, -4.29160589888045, 0.345939603980626,
1.32780224372996, -0.603792378775934, 0.0401273744720227, 0.00122577710465425,
-0.0138137369660334, 1.02828568252612, 0.036, 0, 0, 0, -4.06238321984862,
0.0146353060550562, 0.180980132478463, 1.13664686131079, 0, 0,
0, 0, 1.17627912995327, 0, 0.942517303813378)

init1 <- c(0.00128500849891446, 0.00165353550525411, -0.0457700527122326,
0.852903241638371, 0.0965964850319247, 0.0106681055898985, -0.00833488917465255,
-0.0983272783788632, -3.97083966918167, 0.0224131103967352, 0.143912587571141,
0.829355087608916, -2.06231973874536, -0.197632882886333, -1.70344194695517,
-0.0974725348981226, 0.0782115557466507, -0.00109541725928936,
-0.0654920860767664, 0.755271619906803, 0.036, 0, 0, 0, -4.2102888327682,
-0.00912613468204781, 0.286319221594694, 1.19375723405397, 0,
0, 0, 0, 1.14599024095279, 0, 0.946162518243309)

source('~/Documents/qrmissing/code/BiMLESigma.R')
mod1 <- BiQRGradient(y, R, X, tau = 0.05, niter = 2000,  method = 'heter2', init = init1)
mod3 <- BiQRGradient(y, R, X, tau = 0.3, niter = 2000,  method = 'heter2')
mod5 <- BiQRGradient(y, R, X, tau = 0.5, niter = 2000,  method = 'heter2')
mod7 <- BiQRGradient(y, R, X, tau = 0.7, niter = 2000,  method = 'heter2')
mod9 <- BiQRGradient(y, R, X, tau = 0.95, niter = 2000,  method = 'heter2', init = init9)

coef1 <- coef(mod1)
coef3 <- coef(mod3)
coef5 <- coef(mod5)
coef7 <- coef(mod7)
coef9 <- coef(mod9)

coefw2 <- rbind(coef1[1,], coef3[1, ], coef5[1, ], coef7[1, ], coef9[1,])
coefw3 <- rbind(coef1[2,], coef3[2, ], coef5[2, ], coef7[2, ], coef9[2,])

rownames(coefw2) <- c('tau = 0.05', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.95')
rownames(coefw3) <- c('tau = 0.05', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.95')

colnames(coefw2) <- c('Intercept', 'Age(centered)', 'White', 'BaseWeight')
colnames(coefw3) <- c('Intercept', 'Age(centered)', 'White', 'BaseWeight')

library(xtable)
print(xtable(coefw2))
print(xtable(coefw3))

write.table(rbind(coefw2, coefw3), 'ageracebase.txt', row.names=FALSE)
