## Test BiMLESigma.R
## SIMULATE DATA
rm(list = ls())
source('BiMLESigma.R')
set.seed(100)
p <- 0.5
n <- 400
R <- rbinom(n, 1, p)
b01 <- 1
b00 <- -1
b11 <- 1
b10 <- -1
sigma1 <- 1
sigma0 <- 1

b02 <- 1
b12 <- -1
sigma2 <- 1

g1 <- 0.5
g2 <- 0.3


x <- runif(n, 0, 2)
y <- matrix(0, n, 2)
## for (i in 1:n){
##   if (R[i] == 1){
##     y[i, 1] <- rnorm(1, b01 + x[i] * b11, sigma1)
##     y[i, 2] <- rnorm(1, b02 + x[i] * b12, sigma2)
##   } else {
##     y[i, 1] <- rnorm(1, b00 + x[i] * b10, sigma0)
##     y[i, 2] <- rnorm(1, b00 + x[i] * b10, sigma0)
##   }
## }

for (i in 1:n){
  if (R[i] == 1){
    y[i, 1] <- b01 + x[i] * b11+(1+g1*x[1])*rnorm(1,0, sigma1)
    y[i, 2] <- b02 + x[i] * b12+(1+g2*x[1])*rnorm(1,0, sigma2)
  } else {
    y[i, 1] <- b01 + x[i] * b11+(1+g1*x[1])*rnorm(1,0, sigma0)
    y[i, 2] <- 0
  }
}



## MEAN REGRESSION AND PLOT
## summary(lm(y[,1] ~ x))
## plot(y[,1] ~ x)

## summary(lm(y[,2] ~ x))
## plot(y[,2] ~ x)

## options(warn = 2)

mod5 <- BiQRGradient(y, R, x, tau = 0.5)
print(mod5$param)
plot(y[,1] ~ x)
abline(mod5$param[1:2])

mod7 <- BiQRGradient(y, R, x, tau = 0.7)
abline(mod7$param[1:2])
print(mod7$param)

mod9 <- BiQRGradient(y, R, x, tau = 0.9)
abline(mod9$param[1:2])
print(mod9$param)

mod3 <- BiQRGradient(y, R, x, tau = 0.3)
abline(mod3$param[1:2])
print(mod3$param)

mod1 <- BiQRGradient(y, R, x, tau = 0.1)
abline(mod1$param[1:2])
print(mod1$param)

plot(y[,2][R==1] ~ x[R==1])
abline(mod5$param[7:8])
abline(mod9$param[7:8])
abline(mod7$param[7:8])
abline(mod3$param[7:8])
abline(mod1$param[7:8])


## mod5N <- BiQRGradient(y, R, x, tau = 0.5, niter = 1000, sp = c(1,0,0,1,1))
## print(mod5N$param)

## mod7N <- BiQRGradient(y, R, x, tau = 0.7, niter = 1000, sp = c(1,0,0,1,1))
## print(mod7N$param)

## mod9N <- BiQRGradient(y, R, x, tau = 0.9, niter = 1000, sp = c(1,0,0,1,1))
## print(mod9N$param)

## mod3N <- BiQRGradient(y, R, x, tau = 0.3, niter = 1000, sp = c(1,0,0,1,1))
## print(mod3N$param)

## mod1N <- BiQRGradient(y, R, x, tau = 0.1, niter = 1000, sp = c(1,0,0,1,1))
## print(mod1N$param)

## plot(y[,2][R==1] ~ x[R==1])
## abline(mod5N$param[7:8])
## abline(mod9N$param[7:8])
## abline(mod7N$param[7:8])
## abline(mod3N$param[7:8])
## abline(mod1N$param[7:8])

