## Test BiMLESigma.R
## SIMULATE DATA
rm(list = ls())
source('BiMLESigma.R')
set.seed(1)
p <- 0.5
n <- 200
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

x <- runif(n, 0, 2)
y <- matrix(0, n, 2)
for (i in 1:n){
  if (R[i] == 1){
    y[i, 1] <- rnorm(1, b01 + x[i] * b11, sigma1)
    y[i, 2] <- rnorm(1, b02 + x[i] * b12, sigma2)
  } else {
    y[i, 1] <- rnorm(1, b00 + x[i] * b10, sigma0)
    y[i, 2] <- rnorm(1, b00 + x[i] * b10, sigma0)
  }
}

## MEAN REGRESSION AND PLOT
summary(lm(y[,1] ~ x))
plot(y[,1] ~ x)

summary(lm(y[,2] ~ x))
plot(y[,2] ~ x)

## options(warn = 2)

mod5 <- BiQRGradient(y, R, x, tau = 0.5)
print(mod5$param)
plot(y[,1] ~ x)
abline(mod5$param[1:2])
plot(y[,2][R==1] ~ x[R==1])
abline(mod5$param[7:8])

mod5N <- BiQRGradient(y, R, x, tau = 0.5, niter = 1000, sp = c(1,0,0,1))
print(mod5N$param)
