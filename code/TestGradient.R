## Test Gredient Descent Method
## Time-stamp: <liuminzhao 03/16/2013 20:07:05>

## Test for normal MLE first
## Data
n <- 500
x <- rnorm(n)
beta <- c(1, -1)
sigma <- 2
y <- beta[1] + beta[2] * x + rnorm(n, sd = sigma)

## Target function
LL <- function(beta0, beta1, sigma){
  return(-sum(dnorm(y, mean = beta0 + beta1 * x, sd = sigma, log = T)))
}

cat(LL(beta[1], beta[2], sigma), '\n')

PartialLL0 <- function(beta0, beta1, sigma){
  epsilon <- 0.01
  ll1 <- LL(beta0 + epsilon, beta1, sigma)
  ll0 <- LL(beta0 - epsilon, beta1, sigma)
  return((ll1 - ll0) / 2 / epsilon)
}


PartialLL1 <- function(beta0, beta1, sigma){
  epsilon <- 0.01
  ll1 <- LL(beta0, beta1 + epsilon, sigma)
  ll0 <- LL(beta0, beta1 - epsilon, sigma)
  return((ll1 - ll0) / 2 / epsilon)
}


PartialLLs <- function(beta0, beta1, sigma){
  epsilon <- 0.01
  ll1 <- LL(beta0, beta1, sigma + epsilon)
  ll0 <- LL(beta0, beta1, sigma - epsilon)
  return((ll1 - ll0) / 2 / epsilon)
}


## Begin Gradient Descent
dif <- 1
alpha <- 0.001
beta0 <- beta1 <- sigma <- 4
ll0 <- LL(beta0, beta1, sigma)
llsave <- ll0
iter <- 0
while (dif > 10^-3 & iter < 10000) {
  partial0 <- PartialLL0(beta0, beta1, sigma)
  partial1 <- PartialLL1(beta0, beta1, sigma)
  partials <- PartialLLs(beta0, beta1, sigma)
  beta0 <- beta0 - alpha * partial0
  beta1 <- beta1 - alpha * partial1
  sigma <- sigma - alpha * partials
  ll <- LL(beta0, beta1, sigma)
  dif <- abs(ll - ll0)
  ll0 <- ll
  llsave <- append(llsave, ll)
  iter <- iter + 1
}
cat(beta0, beta1, sigma, '\n')
print(summary(lm(y ~ x)))
l <- length(llsave)
plot(ts(llsave[(l/2): l]))

## pretty successful and quick, need small alpha
