## Test UniMLESigma.R
rm(list = ls())
set.seed(1)
source('UniMLESigma.R')
n <- 500
p <- 0.5
S <- rbinom(n, 1, p)
b01 <- 1
b00 <- -1
b11 <- 1
b10 <- -2
sigma1 <- 1
sigma0 <- 2
x <- runif(n, 0, 2)
y <- rep(0, n)
for (i in 1:n){
  if (S[i] == 1) {
    y[i] <- rnorm(1, b01 + x[i] * b11, sigma1)
  } else {
    y[i] <- rnorm(1, b00 + x[i] * b10, sigma0)
  }
}
y1 <- y[S == 1]
y0 <- y[S == 0]

## mean regression of y on x
summary(lm(y ~ x))
plot(y ~ x)

mod5 <- QRGradient(y, S, x, tau = 0.5)
abline(mod5$param[1:2])
print(mod5$param)

mod7 <- QRGradient(y, S, x, tau = 0.7)
abline(mod7$param[1:2])
print(mod7$param)

mod9 <- QRGradient(y, S, x, tau = 0.9)
abline(mod9$param[1:2])
print(mod9$param)

mod3 <- QRGradient(y, S, x, tau = 0.3)
abline(mod3$param[1:2])
print(mod3$param)

mod1 <- QRGradient(y, S, x, tau = 0.1)
abline(mod1$param[1:2])
print(mod1$param)
