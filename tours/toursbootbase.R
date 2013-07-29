#!/bin/Rscript
## Time-stamp: <liuminzhao 07/29/2013 04:17:14>
## bootstrap on tours data
## weight2 and weight3
## scaled by 1/100
## covariates: 3 treatments and 2 races
## 2013/07/29 add init1 and init9 , modfiied saved file name, change to 10 core
rm(list=ls())
source('~/Documents/qrmissing/code/BiMLESigma.R')
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

  mod1 <- BiQRGradient(y, R, X, tau = 0.05, niter = 500, method = 'heter2', init = init1)
  mod3 <- BiQRGradient(y, R, X, tau = 0.3, niter = 500, method = 'heter2')
  mod5 <- BiQRGradient(y, R, X, tau = 0.5, niter = 500, method = 'heter2')
  mod7 <- BiQRGradient(y, R, X, tau = 0.7, niter = 500, method = 'heter2')
  mod9 <- BiQRGradient(y, R, X, tau = 0.95, niter = 500, method = 'heter2', init = init9)

  coef1 <- coef(mod1)
  coef3 <- coef(mod3)
  coef5 <- coef(mod5)
  coef7 <- coef(mod7)
  coef9 <- coef(mod9)
  count1 <- mod1$converge
  count3 <- mod3$converge
  count5 <- mod5$converge
  count7 <- mod7$converge
  count9 <- mod9$converge

  coefw2 <- c(coef1[1,], coef3[1, ], coef5[1, ], coef7[1, ], coef9[1,])
  coefw3 <- c(coef1[2,], coef3[2, ], coef5[2, ], coef7[2, ], coef9[2,])

  ans <- c(coefw2, coefw3, count1, count3, count5, count7, count9)

}

write.table(result, file = "toursbootagebase.txt", row.names = FALSE, col.names = FALSE)
sendEmail(subject="boot tours base", text="done", address="liuminzhao@gmail.com")

print(proc.time()[3] - start)
