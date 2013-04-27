## Time-stamp: <liuminzhao 04/24/2013 16:50:40>
## bootstrap on tours data
## weight2 and weight3
## scaled by 1/100
## covariates: 3 treatments and 2 races
rm(list=ls())
source('~/Documents/qrmissing/code/BiMLESigma.R')
source('~/Documents/qrmissing/code/sendEmail.R')
library(doMC)
registerDoMC()
options(cores = 8)
set.seed(1)

TOURS <- read.csv('../tours/tours.csv')
TOURS <- subset(TOURS, RACE==1 | RACE==3)
weight2 <- TOURS$wtkg2
weight3 <- TOURS$wtkg3
trt <- TOURS$TREATMENT
race3 <- as.numeric(TOURS$RACE == 3)

## scaled by 1/100
weight2 <- weight2/100
weight3 <- weight3/100

## new dataset
NEWTOURS <- data.frame(w2 = weight2, w3 = weight3, trt, race3)


###############
## BOOTSTRAP 
###############
boot <- 8
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
  X[,2] <- dat$trt == 'O'
  X[,3] <- dat$trt == 'T'
  X[,4] <- dat$race3
  R <- 1 - as.numeric(is.na(dat$w3))
  y[is.na(y[,2]),2] <- 0

  mod1 <- BiQRGradient(y, R, X, tau = 0.1, method = 'heter2')
  mod3 <- BiQRGradient(y, R, X, tau = 0.3, method = 'heter2')
  mod5 <- BiQRGradient(y, R, X, tau = 0.5, method = 'heter2')
  mod7 <- BiQRGradient(y, R, X, tau = 0.7, method = 'heter2')
  mod9 <- BiQRGradient(y, R, X, tau = 0.9, method = 'heter2')

  coef1 <- mysummary(mod1)
  coef3 <- mysummary(mod3)
  coef5 <- mysummary(mod5)
  coef7 <- mysummary(mod7)
  coef9 <- mysummary(mod9)

  coefw2 <- c(coef1[1,], coef3[1, ], coef5[1, ], coef7[1, ], coef9[1,])
  coefw3 <- c(coef1[2,], coef3[2, ], coef5[2, ], coef7[2, ], coef9[2,])

  ans <- c(coefw2, coefw3)
  
}

write.table(result, file = "toursboot.txt", row.names = FALSE, col.names = FALSE)
sendEmail(subject="boot tours", text="done", address="liuminzhao@gmail.com")

print(proc.time()[3] - start)
